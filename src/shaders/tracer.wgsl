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
    roughness: f32,
    metallic: f32,
    transmission: f32,
    ior: f32,
};

struct Triangle {
    v0: vec3f,
    v1: vec3f,
    v2: vec3f,
    color: vec3f,
    emission: vec3f,
    roughness: f32,
    metallic: f32,
    transmission: f32,
    ior: f32,
};

struct HitRecord {
    t: f32,
    p: vec3f,
    normal: vec3f,
    matched: bool,
    color: vec3f,
    emission: vec3f,
    roughness: f32,
    metallic: f32,
    transmission: f32,
    ior: f32,
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

fn random_unit_vector() -> vec3f {
    let z = rand() * 2.0 - 1.0;
    let a = rand() * 2.0 * 3.14159265;
    let r = sqrt(max(0.0, 1.0 - z * z));
    let x = r * cos(a);
    let y = r * sin(a);
    return vec3f(x, y, z);
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
    rec.roughness = s.roughness;
    rec.metallic = s.metallic;
    rec.transmission = s.transmission;
    rec.ior = s.ior;
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
    rec.roughness = tri.roughness;
    rec.metallic = tri.metallic;
    rec.transmission = tri.transmission;
    rec.ior = tri.ior;
    return rec;
}

fn hit_quad(v0: vec3f, v1: vec3f, v2: vec3f, v3: vec3f, color: vec3f, emission: vec3f, roughness: f32, metallic: f32, transmission: f32, ior: f32, r: Ray, t_min: f32, t_max: f32, closest_so_far: f32) -> HitRecord {
    // Quad defined by v0, v1, v2, v3 (CCW or CW).
    // Split into two triangles: v0-v1-v2 and v0-v2-v3
    var final_rec: HitRecord;
    final_rec.matched = false;
    var current_closest = closest_so_far;

    let t1 = Triangle(v0, v1, v2, color, emission, roughness, metallic, transmission, ior);
    let rec1 = hit_triangle(t1, r, t_min, current_closest);
    if (rec1.matched) {
        final_rec = rec1;
        current_closest = rec1.t;
    }

    let t2 = Triangle(v0, v2, v3, color, emission, roughness, metallic, transmission, ior);
    let rec2 = hit_triangle(t2, r, t_min, current_closest);
    if (rec2.matched) {
        final_rec = rec2;
    }
    
    return final_rec;
}


// Scene Helper
fn hit_box(c_min: vec3f, c_max: vec3f, color: vec3f, emission: vec3f, roughness: f32, metallic: f32, transmission: f32, ior: f32, r: Ray, t_min: f32, t_max: f32, closest_so_far: f32) -> HitRecord {
    var rec: HitRecord;
    rec.matched = false;

    let inv_d = 1.0 / r.direction;
    let t0s = (c_min - r.origin) * inv_d;
    let t1s = (c_max - r.origin) * inv_d;

    let t_smaller = min(t0s, t1s);
    let t_larger  = max(t0s, t1s);

    let t_enter = max(max(t_smaller.x, t_smaller.y), t_smaller.z);
    let t_exit  = min(min(t_larger.x, t_larger.y), t_larger.z);

    if (t_exit < t_enter || t_exit < t_min) {
        return rec;
    }
    
    // We want the closest hit in [t_min, closest_so_far]
    // Candidates are t_enter (if > t_min) or t_exit (if inside and t_exit > t_min)
    // Actually standard AABB: if ray origin is outside, hit is t_enter. If inside, hit is t_exit (if we want to see inside faces).
    // Our ray tracer logic usually wants the first hit.
    
    var t_hit = t_enter;
    if (t_hit < t_min) {
        t_hit = t_exit;
        if (t_hit < t_min) {
            return rec;
        }
    }

    if (t_hit >= closest_so_far) {
        return rec;
    }

    // Valid hit
    rec.t = t_hit;
    rec.p = r.origin + r.direction * t_hit;
    rec.matched = true;
    rec.color = color;
    rec.emission = emission;
    rec.roughness = roughness;
    rec.metallic = metallic;
    rec.transmission = transmission;
    rec.ior = ior;

    // Normal calculation
    // A robust way for AABB: compare hit point to bounds
    // Or closer: see which plane t_hit came from.
    // Since we know t_hit equals one of the slab planes...
    let p = rec.p;
    let epsilon = 0.0001;
    var normal = vec3f(0.0);
    
    // Check against faces. Bias slightly outward or use abs diff
    if (abs(p.x - c_min.x) < epsilon) { normal = vec3f(-1.0, 0.0, 0.0); }
    else if (abs(p.x - c_max.x) < epsilon) { normal = vec3f(1.0, 0.0, 0.0); }
    else if (abs(p.y - c_min.y) < epsilon) { normal = vec3f(0.0, -1.0, 0.0); }
    else if (abs(p.y - c_max.y) < epsilon) { normal = vec3f(0.0, 1.0, 0.0); }
    else if (abs(p.z - c_min.z) < epsilon) { normal = vec3f(0.0, 0.0, -1.0); }
    else if (abs(p.z - c_max.z) < epsilon) { normal = vec3f(0.0, 0.0, 1.0); }
    
    // Ensure normal faces against ray? 
    // Usually mathematical normal is outwards.
    // The main loop handles flipping normal if inside (refraction code).
    // So we invoke standard outward normal here.
    rec.normal = normal;

    return rec;
}
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
    
    // Geometry
    
    // Frosted Glass Cuboids & Control (Scaled 0.8x)
    let glass_color = vec3f(1.0, 1.0, 1.0);
    let frost_rough = 0.3;
    let clear_rough = 0.0;
    
    // Y range for all: 0.001 to 0.8 (Lifted slightly to avoid z-fighting)
    // Z range for all: -0.4 to 0.4 (Depth 0.8)
    
    // 1. Frosted (Width 0.32)
    let b1 = hit_box(
        vec3f(-0.68, 0.001, -0.4), vec3f(-0.36, 0.8, 0.4),
        glass_color, vec3f(0.0), frost_rough, 0.0, 1.0, 1.5, r, t_min, t_max, closest_so_far
    );
    if (b1.matched) { closest_so_far = b1.t; final_rec = b1; }

    // 2. Frosted (Width 0.16)
    let b2 = hit_box(
        vec3f(-0.32, 0.001, -0.4), vec3f(-0.16, 0.8, 0.4),
        glass_color, vec3f(0.0), frost_rough, 0.0, 1.0, 1.5, r, t_min, t_max, closest_so_far
    );
    if (b2.matched) { closest_so_far = b2.t; final_rec = b2; }

    // 3. Frosted (Width 0.08)
    let b3 = hit_box(
        vec3f(-0.12, 0.001, -0.4), vec3f(-0.04, 0.8, 0.4),
        glass_color, vec3f(0.0), frost_rough, 0.0, 1.0, 1.5, r, t_min, t_max, closest_so_far
    );
    if (b3.matched) { closest_so_far = b3.t; final_rec = b3; }
    
    // 4. Control Clear (Width 0.8)
    let b4 = hit_box(
        vec3f(0.0, 0.001, -0.4), vec3f(0.8, 0.8, 0.4),
        glass_color, vec3f(0.0), clear_rough, 0.0, 1.0, 1.5, r, t_min, t_max, closest_so_far
    );
    if (b4.matched) { closest_so_far = b4.t; final_rec = b4; }

    // Floor (y=0) - White, Rough
    let rec_floor = hit_quad(
        vec3f(-1.0, 0.0, 1.0), vec3f(1.0, 0.0, 1.0), vec3f(1.0, 0.0, -1.0), vec3f(-1.0, 0.0, -1.0),
        white, vec3f(0.0), 1.0, 0.0, 0.0, 1.5, r, t_min, t_max, closest_so_far
    );
    if (rec_floor.matched) { closest_so_far = rec_floor.t; final_rec = rec_floor; }
    
    // Ceiling (y=2) - White, Rough
    let rec_ceil = hit_quad(
        vec3f(-1.0, 2.0, -1.0), vec3f(1.0, 2.0, -1.0), vec3f(1.0, 2.0, 1.0), vec3f(-1.0, 2.0, 1.0),
        white, vec3f(0.0), 1.0, 0.0, 0.0, 1.5, r, t_min, t_max, closest_so_far
    );
    if (rec_ceil.matched) { closest_so_far = rec_ceil.t; final_rec = rec_ceil; }
    
    // Back Wall (z=-1) - White, Rough
    let rec_back = hit_quad(
        vec3f(-1.0, 0.0, -1.0), vec3f(1.0, 0.0, -1.0), vec3f(1.0, 2.0, -1.0), vec3f(-1.0, 2.0, -1.0),
        white, vec3f(0.0), 1.0, 0.0, 0.0, 1.5, r, t_min, t_max, closest_so_far
    );
    if (rec_back.matched) { closest_so_far = rec_back.t; final_rec = rec_back; }
    
    // Left Wall (x=-1) - Red, Rough
    let rec_left = hit_quad(
        vec3f(-1.0, 0.0, 1.0), vec3f(-1.0, 0.0, -1.0), vec3f(-1.0, 2.0, -1.0), vec3f(-1.0, 2.0, 1.0),
        red, vec3f(0.0), 1.0, 0.0, 0.0, 1.5, r, t_min, t_max, closest_so_far
    );
    if (rec_left.matched) { closest_so_far = rec_left.t; final_rec = rec_left; }
    
    // Right Wall (x=1) - Green, Rough
    let rec_right = hit_quad(
        vec3f(1.0, 0.0, -1.0), vec3f(1.0, 0.0, 1.0), vec3f(1.0, 2.0, 1.0), vec3f(1.0, 2.0, -1.0),
        green, vec3f(0.0), 1.0, 0.0, 0.0, 1.5, r, t_min, t_max, closest_so_far
    );
    if (rec_right.matched) { closest_so_far = rec_right.t; final_rec = rec_right; }
    
    // Light (Ceiling patch)
    let light_size = 0.5;
    let rec_light = hit_quad(
        vec3f(-light_size/2.0, 1.99, -light_size/2.0), vec3f(light_size/2.0, 1.99, -light_size/2.0),
        vec3f(light_size/2.0, 1.99, light_size/2.0), vec3f(-light_size/2.0, 1.99, light_size/2.0),
        vec3f(0.0), light, 1.0, 0.0, 0.0, 1.5, r, t_min, t_max, closest_so_far
    );
    if (rec_light.matched) { closest_so_far = rec_light.t; final_rec = rec_light; }


    
    return final_rec;
}

// PBR Helpers
fn ortho_basis(n: vec3f) -> mat3x3f {
    let f = normalize(n);
    var r = cross(vec3f(0.0, 1.0, 0.0), f);
    if (length(r) < 0.1) {
        r = cross(vec3f(1.0, 0.0, 0.0), f);
    }
    r = normalize(r);
    let u = cross(f, r);
    return mat3x3f(r, u, f);
}

fn schlick_fresnel(u: f32, f0: vec3f) -> vec3f {
    return f0 + (vec3f(1.0) - f0) * pow(1.0 - u, 5.0);
}

fn sample_ggx(n: vec3f, roughness: f32, r1: f32, r2: f32) -> vec3f {
    let a = roughness * roughness;
    let phi = 2.0 * 3.14159 * r1;
    let cos_theta = sqrt((1.0 - r2) / (1.0 + (a * a - 1.0) * r2));
    let sin_theta = sqrt(1.0 - cos_theta * cos_theta);

    let h = vec3f(sin_theta * cos(phi), sin_theta * sin(phi), cos_theta); // Tangent space
    
    let basis = ortho_basis(n);
    return normalize(basis * h); // World space
}

fn refract(uv: vec3f, n: vec3f, etai_over_etat: f32) -> vec3f {
    let cos_theta = min(dot(-uv, n), 1.0);
    let r_out_perp = etai_over_etat * (uv + cos_theta * n);
    let r_out_parallel = -sqrt(abs(1.0 - dot(r_out_perp, r_out_perp))) * n;
    return r_out_perp + r_out_parallel;
}

fn ray_color(r_in: Ray) -> vec3f {
    var r = r_in;
    var attenuation = vec3f(1.0);
    var accumulated_light = vec3f(0.0);

    // Bounce limit
    for (var i = 0u; i < 10u; i++) {
        let rec = hit_world(r, 0.001, 10000.0);
        if (rec.matched) {
            // Emitted light
            accumulated_light += attenuation * rec.emission;
            
            // Material properties
            let normal = rec.normal;
            let view_dir = normalize(-r.direction);
            var albedo = rec.color;
            let roughness = rec.roughness;
            let metallic = rec.metallic;
            let transmission = rec.transmission;
            let ior = rec.ior;
            
            // Variables for scatter calc
            var scatter_dir: vec3f;
            var current_attenuation = vec3f(1.0);
            
            // Detect if we are inside the object
            let front_face = dot(r.direction, rec.normal) < 0.0;
            var real_normal = rec.normal;
            var refraction_ratio = 1.0 / ior;
            if (!front_face) {
                real_normal = -rec.normal;
                refraction_ratio = ior;
            }

            // Probability for lobes
            // Simple priority: Metallic -> Transmissive -> Diffuse
            
            // Fresnel (Schlick)
            // Use base F0 for dielectric (0.04) or albedo (for metal)
            var f0 = vec3f(0.04); 
            f0 = mix(f0, albedo, metallic);
            
            let cos_theta = min(dot(-r.direction, real_normal), 1.0);
            let F = schlick_fresnel(cos_theta, f0);
            
            // Reflection Probability via Fresnel
            // For dielectrics, F determines Reflect vs Refract
            let p_reflect_fresnel = (F.x + F.y + F.z) / 3.0;
            
            // 1. Check for Metal
            if (metallic > 0.5) {
                 // Metal is always specular reflection
                 // Roughness jitter
                 let h = sample_ggx(real_normal, roughness, rand(), rand());
                 scatter_dir = reflect(-view_dir, h);
                 current_attenuation = albedo; // Metals absorb based on color
            } 
            // 2. Check for Transmission (Glass)
            else if (transmission > 0.0) {
                // Glass / Dielectric
                
                // 1. Sample Microfacet Normal (h)
                let h = sample_ggx(real_normal, roughness, rand(), rand());
                
                // 2. Calculate Fresnel using h and view direction
                // Ensure dot product is positive (view and h on same side conceptually for F)
                let v_dot_h = dot(view_dir, h);
                let F = schlick_fresnel(abs(v_dot_h), f0);
                let p_reflect = (F.x + F.y + F.z) / 3.0;
                
                // 3. Russian Roulette: Reflect or Refract
                let do_reflect = (rand() < p_reflect);
                
                // Beer's Law: Apply absorption if we just traveled THROUGH the medium
                // We are 'inside' if we are hitting the back face (!front_face)
                var absorption_factor = vec3f(1.0);
                if (!front_face) {
                    // Simple Beer's Law approximation
                    // Absorbance derived from albedo: A = -log(albedo)
                    // Density factor to tune the "filled" look
                    let density = 2.0;
                    let absorbance = -log(max(albedo, vec3f(0.001))) * density;
                    absorption_factor = exp(-absorbance * rec.t);
                }
                
                if (do_reflect) {
                     // Microfacet Reflection
                     scatter_dir = reflect(-view_dir, h);
                     current_attenuation = absorption_factor; // Apply absorption from path
                } else {
                    // Refract
                    // Try to refract using microfacet normal h
                    let refraction_dir = refract(-view_dir, h, refraction_ratio);
                    
                    if (length(refraction_dir) == 0.0) {
                        // Total Internal Reflection (TIR) -> Must Reflect
                        scatter_dir = reflect(-view_dir, h);
                        current_attenuation = absorption_factor;
                    } else {
                        scatter_dir = normalize(refraction_dir);
                        // On Transmission, do we apply albedo tint at surface?
                        // Physically, absorption happens in volume.
                        // Standard practice: pure surface transmission is 1.0 (minus fresnel which we handled)
                        // But we apply volume absorption.
                        current_attenuation = absorption_factor;
                    }
                }
            }
            // 3. Diffuse / Glossy Plastic
            else {
                // Dielectric (Plastic/Matte)
                // Mix Specular and Diffuse based on Fresnel
                let prob_spec = p_reflect_fresnel; 
                let sample_specular = (rand() < prob_spec);
                
                if (sample_specular) {
                     let h = sample_ggx(real_normal, roughness, rand(), rand());
                     scatter_dir = reflect(-view_dir, h);
                     current_attenuation = vec3f(1.0); // White highlight
                } else {
                    // Diffuse
                    let diffuse_target = rec.p + real_normal + random_unit_vector();
                    scatter_dir = normalize(diffuse_target - rec.p);
                    current_attenuation = albedo;
                }
            }
            
            r.origin = rec.p;
            r.direction = scatter_dir;
            attenuation *= current_attenuation;
            
            // Russian Roulette for termination
            if (length(attenuation) < 0.001) {
                break;
            }
        } else {
            // Sky
            accumulated_light += vec3f(0.0);
            break; 
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
