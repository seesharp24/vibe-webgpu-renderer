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
    materialType: u32, // 0: Lambertian, 1: Metal, 2: Dielectric
};

struct Triangle {
    v0: vec3f,
    v1: vec3f,
    v2: vec3f,
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
    raysPerFrame: u32,
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

fn hit_triangle(tri: Triangle, r: Ray, t_min: f32, t_max: f32) -> HitRecord {
    var rec: HitRecord;
    rec.matched = false;

    let v0v1 = tri.v1 - tri.v0;
    let v0v2 = tri.v2 - tri.v0;
    let pvec = cross(r.direction, v0v2);
    let det = dot(v0v1, pvec);

    // Backface culling if desired, but for cornell box walls we want two-sided or correct normals
    // if det is near zero, parallel
    if (abs(det) < 0.00001) { return rec; }

    let inv_det = 1.0 / det;
    let tvec = r.origin - tri.v0;
    let u = dot(tvec, pvec) * inv_det;

    if (u < 0.0 || u > 1.0) { return rec; }

    let qvec = cross(tvec, v0v1);
    let v = dot(r.direction, qvec) * inv_det;

    if (v < 0.0 || u + v > 1.0) { return rec; }

    let t = dot(v0v2, qvec) * inv_det;

    if (t < t_min || t >= t_max) { return rec; }

    rec.t = t;
    rec.p = r.origin + r.direction * rec.t;
    // Normal:
    // For single sided: cross(v0v1, v0v2)
    // We want normal to face the ray
    var normal = normalize(cross(v0v1, v0v2));
    if (dot(normal, r.direction) > 0.0) {
        normal = -normal;
    }
    rec.normal = normal;
    
    rec.matched = true;
    rec.color = tri.color;
    rec.emission = tri.emission;
    return rec;
}

fn hit_quad(v0: vec3f, v1: vec3f, v2: vec3f, v3: vec3f, color: vec3f, emission: vec3f, r: Ray, t_min: f32, t_max: f32, closest_so_far: f32) -> HitRecord {
    // Quad defined by v0, v1, v2, v3 (CCW or CW).
    // Split into two triangles: v0-v1-v2 and v0-v2-v3
    var final_rec: HitRecord;
    final_rec.matched = false;
    var current_closest = closest_so_far;

    let t1 = Triangle(v0, v1, v2, color, emission);
    let rec1 = hit_triangle(t1, r, t_min, current_closest);
    if (rec1.matched) {
        final_rec = rec1;
        current_closest = rec1.t;
    }

    let t2 = Triangle(v0, v2, v3, color, emission);
    let rec2 = hit_triangle(t2, r, t_min, current_closest);
    if (rec2.matched) {
        final_rec = rec2;
    }
    
    return final_rec;
}


// Scene Helper
fn hit_world(r: Ray, t_min: f32, t_max: f32) -> HitRecord {
    var hit_anything = false;
    var closest_so_far = t_max;
    var final_rec: HitRecord;
    final_rec.matched = false;

    // Hardcoded Scene
    // Cornell Box Scene
    // Scale: -1 to 1 in X, 0 to 2 in Y, -1 to 1 in Z
    
    // Materials
    let red = vec3f(0.65, 0.05, 0.05);
    let white = vec3f(0.73, 0.73, 0.73);
    let green = vec3f(0.12, 0.45, 0.15);
    let light = vec3f(15.0, 15.0, 15.0);
    
    // Geometry
    
    // Floor (y=0)
    // v0(-1,0,1), v1(1,0,1), v2(1,0,-1), v3(-1,0,-1)
    let rec_floor = hit_quad(
        vec3f(-1.0, 0.0, 1.0), vec3f(1.0, 0.0, 1.0), vec3f(1.0, 0.0, -1.0), vec3f(-1.0, 0.0, -1.0),
        white, vec3f(0.0), r, t_min, t_max, closest_so_far
    );
    if (rec_floor.matched) { closest_so_far = rec_floor.t; final_rec = rec_floor; }
    
    // Ceiling (y=2)
    let rec_ceil = hit_quad(
        vec3f(-1.0, 2.0, -1.0), vec3f(1.0, 2.0, -1.0), vec3f(1.0, 2.0, 1.0), vec3f(-1.0, 2.0, 1.0),
        white, vec3f(0.0), r, t_min, t_max, closest_so_far
    );
    if (rec_ceil.matched) { closest_so_far = rec_ceil.t; final_rec = rec_ceil; }
    
    // Back Wall (z=-1)
    let rec_back = hit_quad(
        vec3f(-1.0, 0.0, -1.0), vec3f(1.0, 0.0, -1.0), vec3f(1.0, 2.0, -1.0), vec3f(-1.0, 2.0, -1.0),
        white, vec3f(0.0), r, t_min, t_max, closest_so_far
    );
    if (rec_back.matched) { closest_so_far = rec_back.t; final_rec = rec_back; }
    
    // Left Wall (x=-1) - Red
    let rec_left = hit_quad(
        vec3f(-1.0, 0.0, 1.0), vec3f(-1.0, 0.0, -1.0), vec3f(-1.0, 2.0, -1.0), vec3f(-1.0, 2.0, 1.0),
        red, vec3f(0.0), r, t_min, t_max, closest_so_far
    );
    if (rec_left.matched) { closest_so_far = rec_left.t; final_rec = rec_left; }
    
    // Right Wall (x=1) - Green
    let rec_right = hit_quad(
        vec3f(1.0, 0.0, -1.0), vec3f(1.0, 0.0, 1.0), vec3f(1.0, 2.0, 1.0), vec3f(1.0, 2.0, -1.0),
        green, vec3f(0.0), r, t_min, t_max, closest_so_far
    );
    if (rec_right.matched) { closest_so_far = rec_right.t; final_rec = rec_right; }
    
    // Light (Ceiling patch)
    let light_size = 0.5;
    let rec_light = hit_quad(
        vec3f(-light_size/2.0, 1.99, -light_size/2.0), vec3f(light_size/2.0, 1.99, -light_size/2.0),
        vec3f(light_size/2.0, 1.99, light_size/2.0), vec3f(-light_size/2.0, 1.99, light_size/2.0),
        vec3f(0.0), light, r, t_min, t_max, closest_so_far
    );
    if (rec_light.matched) { closest_so_far = rec_light.t; final_rec = rec_light; }

    // Sphere 1 (Mirror-like)
    let s1 = Sphere(vec3f(-0.4, 0.4, -0.3), 0.4, vec3f(0.9, 0.9, 0.9), vec3f(0.0), 0u); 
    // Just putting generic spheres back for now, we don't have materials implemented well yet in hit_world structs fully
    // but the Sphere struct was updated earlier to have materialType, but the Sphere struct in hit_world was using the old def.
    // Wait, the Sphere definition in line 111 needs to be updated.
    
    let rec_s1 = hit_sphere(s1, r, t_min, closest_so_far);
    if (rec_s1.matched) { closest_so_far = rec_s1.t; final_rec = rec_s1; }

    // Sphere 2
    let s2 = Sphere(vec3f(0.4, 0.4, 0.3), 0.4, vec3f(0.2, 0.2, 0.8), vec3f(0.0), 0u);
    let rec_s2 = hit_sphere(s2, r, t_min, closest_so_far);
    if (rec_s2.matched) { closest_so_far = rec_s2.t; final_rec = rec_s2; }
    
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
    
    let r_base = Ray(ray_origin, ray_direction);
    
    // Multi-sample loop
    var accumulated_color = vec3f(0.0);
    for (var i = 0u; i < scene.raysPerFrame; i++) {
        // Init RNG for each sub-sample
        init_rng(global_id.xy, scene.frameCount + i); // Offset seed effectively
        
        // Jitter for anti-aliasing (optional, but good for convergence)
        // Re-calculating uv with jitter would be better, but for now let's just trace
        // standard path with different random numbers in ray_color
        
        // Note: ray_color uses global state 'seed' which is modified by rand() calls.
        // We need to re-seed or let it continue? 
        // Best approach: Just call ray_color multiple times.
        // But ray_color modifies state. So subsequent calls are random.
        // However, 'u' and 'v' above are constant for this pixel.
        // Ideally we should jitter u/v inside the loop for AA.
        // For simpler implementation now: just jitter the path (ray_color does that).
        
        accumulated_color += ray_color(r_base);
    }
    
    var pixel_color = accumulated_color / f32(scene.raysPerFrame);

    // Accumulation
    // Read from historyTex (texture_2d, need load with mip level 0, coords are integers)
    let old_color = textureLoad(historyTex, global_id.xy, 0);
    let frame_count = f32(scene.frameCount);
    
    var final_color = pixel_color;
    if (frame_count > 0.0) {
       // frame_count is previous total rays.
       // We just added raysPerFrame rays.
       // Current total = frame_count + raysPerFrame
       // Weight for new batch = raysPerFrame / (frame_count + raysPerFrame)
       let weight = f32(scene.raysPerFrame) / (frame_count + f32(scene.raysPerFrame));
       final_color = mix(old_color.rgb, pixel_color, weight);
    }
    
    // Write to outputTex (storage)
    textureStore(outputTex, global_id.xy, vec4f(final_color, 1.0));
}
