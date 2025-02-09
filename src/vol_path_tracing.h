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
Spectrum vol_path_tracing_2([[maybe_unused]]const Scene &scene,
                            [[maybe_unused]]int x, [[maybe_unused]]int y, /* pixel coordinates */
                            [[maybe_unused]]pcg32_state &rng) {
    // Homework 2: implememt this!
    return make_zero_spectrum();
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
