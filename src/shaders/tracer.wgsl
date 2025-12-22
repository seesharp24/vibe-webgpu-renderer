// tracer.wgsl

struct Ray {
    origin: vec3f,
    direction: vec3f,
};

struct Sphere {
    center: vec3f,
    radius: f32,
    color: vec3f,
    emission: vec3f,
};

struct HitRecord {
    t: f32,
    p: vec3f,
    normal: vec3f,
    matched: bool,
    color: vec3f,
    emission: vec3f,
};

struct SceneUniforms {
    cameraPosition: vec3f,
    cameraForward: vec3f,
    cameraRight: vec3f,
    cameraUp: vec3f,
    frameCount: u32,
    resolution: vec2f,
};

@group(0) @binding(0) var<uniform> scene: SceneUniforms;
@group(0) @binding(1) var historyTex: texture_2d<f32>;
@group(0) @binding(2) var outputTex: texture_storage_2d<rgba32float, write>;

// RNG State
var<private> seed: u32;

fn init_rng(pixel: vec2u, frame: u32) {
    seed = (pixel.x * 1973u + pixel.y * 9277u + frame * 26699u) | 1u;
}

fn rand() -> f32 {
    seed = (seed ^ 61u) ^ (seed >> 16u);
    seed = seed * 9u;
    seed = seed ^ (seed >> 4u);
    seed = seed * 668265261u;
    seed = seed ^ (seed >> 15u);
    return f32(seed) / 4294967296.0;
}

fn random_in_unit_sphere() -> vec3f {
    for (var i = 0u; i < 10u; i++) {
        let p = vec3f(rand(), rand(), rand()) * 2.0 - vec3f(1.0);
        if (dot(p, p) < 1.0) {
            return p;
        }
    }
    return normalize(vec3f(rand() - 0.5, rand() - 0.5, rand() - 0.5));
}

fn random_unit_vector() -> vec3f {
    return normalize(random_in_unit_sphere());
}

// Intersections

fn hit_sphere(s: Sphere, r: Ray, t_min: f32, t_max: f32) -> HitRecord {
    var rec: HitRecord;
    rec.matched = false;

    let oc = r.origin - s.center;
    let a = dot(r.direction, r.direction);
    let half_b = dot(oc, r.direction);
    let c = dot(oc, oc) - s.radius * s.radius;
    let discriminant = half_b * half_b - a * c;

    if (discriminant < 0.0) {
        return rec;
    }
    let sqrtd = sqrt(discriminant);

    // Find nearest root in range
    var root = (-half_b - sqrtd) / a;
    if (root <= t_min || root >= t_max) {
        root = (-half_b + sqrtd) / a;
        if (root <= t_min || root >= t_max) {
            return rec;
        }
    }

    rec.t = root;
    rec.p = r.origin + r.direction * rec.t;
    rec.normal = (rec.p - s.center) / s.radius;
    rec.matched = true;
    rec.color = s.color;
    rec.emission = s.emission;
    return rec;
}

// Scene Helper
fn hit_world(r: Ray, t_min: f32, t_max: f32) -> HitRecord {
    var hit_anything = false;
    var closest_so_far = t_max;
    var final_rec: HitRecord;
    final_rec.matched = false;

    // Hardcoded Scene
    // 1. Center Sphere
    let s1 = Sphere(vec3f(0.0, 0.0, 0.0), 1.0, vec3f(0.5, 0.1, 0.1), vec3f(0.0)); // Reddish
    let rec1 = hit_sphere(s1, r, t_min, closest_so_far);
    if (rec1.matched) {
        hit_anything = true;
        closest_so_far = rec1.t;
        final_rec = rec1;
    }

    // 2. Ground Sphere (huge)
    let s2 = Sphere(vec3f(0.0, -100.5, 0.0), 100.0, vec3f(0.5, 0.8, 0.5), vec3f(0.0)); // Greenish ground
    let rec2 = hit_sphere(s2, r, t_min, closest_so_far);
    if (rec2.matched) {
        hit_anything = true;
        closest_so_far = rec2.t;
        final_rec = rec2;
    }
    
    // 3. Light Sphere
    let s3 = Sphere(vec3f(0.0, 1.5, 0.0), 0.5, vec3f(0.0, 0.0, 0.0), vec3f(10.0, 10.0, 10.0)); // Bright white light
    let rec3 = hit_sphere(s3, r, t_min, closest_so_far);
    if (rec3.matched) {
        hit_anything = true;
        closest_so_far = rec3.t;
        final_rec = rec3;
    }
    
    return final_rec;
}

fn ray_color(r_in: Ray) -> vec3f {
    var r = r_in;
    var attenuation = vec3f(1.0);
    var accumulated_light = vec3f(0.0);

    // Bounce limit
    for (var i = 0u; i < 5u; i++) {
        let rec = hit_world(r, 0.001, 10000.0);
        if (rec.matched) {
            // Emitted light from hit object
            accumulated_light += attenuation * rec.emission;
            
            // Lambertian Scatter
            let scatter_target = rec.p + rec.normal + random_unit_vector();
            
            var direction = scatter_target - rec.p;
            if (length(direction) < 0.001) {
                direction = rec.normal;
            }
            
            r.origin = rec.p;
            r.direction = normalize(direction);
            
            attenuation *= rec.color;
            
            // Russian Roulette
            if (length(attenuation) < 0.001) {
                break;
            }
        } else {
            // Sky background (Daylight)
            let unit_direction = normalize(r.direction);
            let t = 0.5 * (unit_direction.y + 1.0);
            let sky = ((1.0 - t) * vec3f(1.0, 1.0, 1.0) + t * vec3f(0.5, 0.7, 1.0)) * 0.05;
            
            accumulated_light += attenuation * sky;
            break; // Exit loop
        }
    }
    return accumulated_light;
}

@compute @workgroup_size(8, 8)
fn main(@builtin(global_invocation_id) global_id: vec3u) {
    let resolution = vec2u(scene.resolution);
    if (global_id.x >= resolution.x || global_id.y >= resolution.y) {
        return;
    }

    init_rng(global_id.xy, scene.frameCount);

    // Standard camera ray gen
    let aspect_ratio = f32(resolution.x) / f32(resolution.y);
    let u = (f32(global_id.x) + rand()) / f32(resolution.x);
    let v = (f32(resolution.y - 1u - global_id.y) + rand()) / f32(resolution.y); // Flip Y
    
    let viewport_height = 2.0;
    let viewport_width = viewport_height * aspect_ratio;
    
    let horizontal = scene.cameraRight * viewport_width;
    let vertical = scene.cameraUp * viewport_height;
    let lower_left_corner = scene.cameraPosition - horizontal/2.0 - vertical/2.0 + scene.cameraForward;
    
    let ray_origin = scene.cameraPosition;
    let ray_direction = normalize(lower_left_corner + u*horizontal + v*vertical - ray_origin);
    
    let r = Ray(ray_origin, ray_direction);
    var pixel_color = ray_color(r);

    // Accumulation
    // Read from historyTex (texture_2d, need load with mip level 0, coords are integers)
    let old_color = textureLoad(historyTex, global_id.xy, 0);
    let frame_count = f32(scene.frameCount);
    
    var final_color = pixel_color;
    if (frame_count > 0.0) {
       let weight = 1.0 / (frame_count + 1.0);
       final_color = mix(old_color.rgb, pixel_color, weight);
    }
    
    // Write to outputTex (storage)
    textureStore(outputTex, global_id.xy, vec4f(final_color, 1.0));
}
