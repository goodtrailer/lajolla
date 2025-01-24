#include "../microfacet.h"

Spectrum eval_op::operator()(const DisneySheen &bsdf) const {
    if (dot(vertex.geometric_normal, dir_in) < 0 ||
            dot(vertex.geometric_normal, dir_out) < 0) {
        // No light below the surface
        return make_zero_spectrum();
    }
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) < 0) {
        frame = -frame;
    }

    Spectrum base_color = eval(bsdf.base_color, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real sheen_tint = eval(bsdf.sheen_tint, vertex.uv, vertex.uv_screen_size, texture_pool);

    Vector3 dir_half = normalize(dir_in + dir_out);

    Real dot_ho = abs(dot(dir_half, dir_out));
    Real dot_no = abs(dot(frame.n, dir_out));
    Real comp_ho = 1 - dot_ho;

    Real l = luminance(base_color);
    Spectrum C_tint = l > 0 ? base_color / l : make_const_spectrum(1);
    Spectrum C_sheen = (1 - sheen_tint) + sheen_tint * C_tint;

    Spectrum f = C_sheen * comp_ho * comp_ho * comp_ho * comp_ho * comp_ho * dot_no;

    return f;
}

Real pdf_sample_bsdf_op::operator()(const DisneySheen &) const {
    if (dot(vertex.geometric_normal, dir_in) < 0 ||
            dot(vertex.geometric_normal, dir_out) < 0) {
        // No light below the surface
        return 0;
    }
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) < 0) {
        frame = -frame;
    }

    // Same as Lambertian, we importance sample the cosine hemisphere domain.
    // Doesn't 100% match the BSDF, but it's pretty close and analytically simple.
    return fmax(dot(frame.n, dir_out), Real(0)) / c_PI;
}

std::optional<BSDFSampleRecord>
        sample_bsdf_op::operator()(const DisneySheen &) const {
    if (dot(vertex.geometric_normal, dir_in) < 0) {
        // No light below the surface
        return {};
    }
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) < 0) {
        frame = -frame;
    }

    Vector3 dir_out = to_world(frame, sample_cos_hemisphere(rnd_param_uv));

    return BSDFSampleRecord { dir_out, 0, 1 };
}

TextureSpectrum get_texture_op::operator()(const DisneySheen &bsdf) const {
    return bsdf.base_color;
}
