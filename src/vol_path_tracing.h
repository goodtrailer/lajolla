#pragma once

// The simplest volumetric renderer: 
// single absorption only homogeneous volume
// only handle directly visible light sources
Spectrum vol_path_tracing_1(const Scene &scene,
                            int x, int y, /* pixel coordinates */
                            pcg32_state &rng) {
    constexpr RayDifferential ray_diff { Real(0), Real(0) };

    int w = scene.camera.width;
    int h = scene.camera.height;
    Vector2 screen_pos((x + next_pcg32_real<Real>(rng)) / w,
        (y + next_pcg32_real<Real>(rng)) / h);
    Ray ray = sample_primary(scene.camera, screen_pos);

    Spectrum radiance = make_zero_spectrum();

    if (auto vertex = intersect(scene, ray, ray_diff)) {
        if (is_light(scene.shapes[vertex->shape_id])) {
            Real sigma_a = get_sigma_a(scene.media[vertex->exterior_medium_id], vertex->position).x;
            Real t = length(vertex->position - ray.org);

            Real transmission = exp(-sigma_a * t);
            Spectrum Le = emission(*vertex, -ray.dir, scene);
            radiance = transmission * Le;
        }
    }

    return radiance;
}

// The second simplest volumetric renderer: 
// single monochromatic homogeneous volume with single scattering,
// no need to handle surface lighting, only directly visible light source
Spectrum vol_path_tracing_2(const Scene &scene,
                            int x, int y, /* pixel coordinates */
                            pcg32_state &rng) {
    constexpr RayDifferential ray_diff { Real(0), Real(0) };

    int w = scene.camera.width;
    int h = scene.camera.height;
    Vector2 screen_pos((x + next_pcg32_real<Real>(rng)) / w,
        (y + next_pcg32_real<Real>(rng)) / h);
    Ray ray = sample_primary(scene.camera, screen_pos);
    
    const Medium& medium = scene.media[scene.camera.medium_id];

    Real sigma_a = get_sigma_a(medium, ray.org).x;
    Real sigma_s = get_sigma_s(medium, ray.org).x;
    Real sigma_t = sigma_a + sigma_s;

    Real u = next_pcg32_real<Real>(rng);
    Real t = -log(1 - u) / sigma_t;
    
    std::optional<PathVertex> hit = intersect(scene, ray, ray_diff);
    Real t_hit = hit ? length(hit->position - ray.org) : infinity<Real>();

    Spectrum radiance = make_zero_spectrum();

    if (t < t_hit) {
        Real pdf_t = exp(-sigma_t * t) * sigma_t;
        Real transmittance = exp(-sigma_t * t);

        Vector3 p = ray.org + t * ray.dir;

        Spectrum L_scatter1 = make_zero_spectrum();
        {
            int light_id = sample_light(scene, next_pcg32_real<Real>(rng));
            const Light& light = scene.lights[light_id];
            Vector2 light_uv { next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng) };
            Real light_w = next_pcg32_real<Real>(rng);

            auto [p_, n_] = sample_point_on_light(light, p, light_uv, light_w, scene);
            Vector3 dir_light = normalize(p_ - p);

            Real pdf_p_ = light_pmf(scene, light_id) * pdf_point_on_light(light, { p_, n_ }, p, scene);
            Real transmittance_ = exp(-sigma_t * length(p - p_));
            
            Spectrum phase = eval(get_phase_function(medium), -ray.dir, dir_light);
            Spectrum L_e = emission(light, -dir_light, 0, { p_, n_ }, scene);

            Real dwdp_ = 0;

            Real epsilon = get_shadow_epsilon(scene);
            Ray shadow_ray { p, dir_light, epsilon, (1 - epsilon) * length(p - p_) };
            if (!occluded(scene, shadow_ray)) {
                dwdp_ = abs(dot(dir_light, n_)) / length_squared(p - p_);
            }
            
            L_scatter1 = phase * L_e * transmittance_ * dwdp_ / pdf_p_;
        }

        radiance = sigma_s * L_scatter1 * transmittance / pdf_t;
    } else {
        Real pr_t = exp(-sigma_t * t_hit);
        Real transmittance = exp(-sigma_t * t_hit);

        Spectrum L_e = make_zero_spectrum();
        if (hit && is_light(scene.shapes[hit->shape_id])) {
            L_e = emission(*hit, -ray.dir, scene);
        }

        radiance = L_e * transmittance / pr_t;
    }

    return radiance;
}

// The third volumetric renderer (not so simple anymore): 
// multiple monochromatic homogeneous volumes with multiple scattering
// no need to handle surface lighting, only directly visible light source
Spectrum vol_path_tracing_3([[maybe_unused]]const Scene &scene,
                            [[maybe_unused]]int x, [[maybe_unused]]int y, /* pixel coordinates */
                            [[maybe_unused]]pcg32_state &rng) {
    // Homework 2: implememt this!
    return make_zero_spectrum();
}

// The fourth volumetric renderer: 
// multiple monochromatic homogeneous volumes with multiple scattering
// with MIS between next event estimation and phase function sampling
// still no surface lighting
Spectrum vol_path_tracing_4([[maybe_unused]]const Scene &scene,
                            [[maybe_unused]]int x, [[maybe_unused]]int y, /* pixel coordinates */
                            [[maybe_unused]]pcg32_state &rng) {
    // Homework 2: implememt this!
    return make_zero_spectrum();
}

// The fifth volumetric renderer: 
// multiple monochromatic homogeneous volumes with multiple scattering
// with MIS between next event estimation and phase function sampling
// with surface lighting
Spectrum vol_path_tracing_5([[maybe_unused]]const Scene &scene,
                            [[maybe_unused]]int x, [[maybe_unused]]int y, /* pixel coordinates */
                            [[maybe_unused]]pcg32_state &rng) {
    // Homework 2: implememt this!
    return make_zero_spectrum();
}

// The final volumetric renderer: 
// multiple chromatic heterogeneous volumes with multiple scattering
// with MIS between next event estimation and phase function sampling
// with surface lighting
Spectrum vol_path_tracing([[maybe_unused]]const Scene &scene,
                          [[maybe_unused]]int x, [[maybe_unused]]int y, /* pixel coordinates */
                          [[maybe_unused]]pcg32_state &rng) {
    // Homework 2: implememt this!
    return make_zero_spectrum();
}
