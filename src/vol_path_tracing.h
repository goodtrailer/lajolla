#pragma once

// The simplest volumetric renderer: 
// single absorption only homogeneous volume
// only handle directly visible light sources
Spectrum vol_path_tracing_1([[maybe_unused]]const Scene &scene,
                            [[maybe_unused]]int x, [[maybe_unused]]int y, /* pixel coordinates */
                            [[maybe_unused]]pcg32_state &rng) {
    // Homework 2: implememt this!
    return make_zero_spectrum();
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
