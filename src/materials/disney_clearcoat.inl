#include "../microfacet.h"

Spectrum eval_op::operator()(const DisneyClearcoat &bsdf) const {
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) < 0) {
        frame = -frame;
    }

    Real clearcoat_gloss = eval(bsdf.clearcoat_gloss, vertex.uv, vertex.uv_screen_size, texture_pool);

    Vector3 dir_half = normalize(dir_in + dir_out);

    Real dot_ni = abs(dot(frame.n, dir_in));
    Real dot_nh = abs(dot(frame.n, dir_half));
    Real dot_ho = abs(dot(dir_half, dir_out));
    Real comp_ho = 1 - dot_ho;

    Real fresnel;
    {
        constexpr Real ior = 1.5;
        constexpr Real ior_air = 1;
        constexpr Real fresnel_90 = (ior - ior_air) * (ior - ior_air) / ((ior + ior_air) * (ior + ior_air));

        fresnel = fresnel_90 + (1 - fresnel_90) * comp_ho * comp_ho * comp_ho * comp_ho * comp_ho;
    }

    Real distribution;
    {
        Real alpha = (1 - clearcoat_gloss) * Real(0.1) + clearcoat_gloss * Real(0.001);
        Real alpha2 = alpha * alpha;

        distribution = (alpha2 - 1) / (c_PI * log(alpha2) * (1 + (alpha2 - 1) * dot_nh * dot_nh));
    }

    Real masking;
    {
        constexpr Real roughness = Real(0.5);
        constexpr Real alpha = roughness * roughness;

        Vector3 dir_local_in = to_local(frame, dir_in);
        Vector3 dir_local_out = to_local(frame, dir_out);

        Vector3 aux_in = Vector3 { alpha, alpha, Real(1) } * dir_local_in;
        Vector3 aux_out = Vector3 { alpha, alpha, Real(1) } * dir_local_out;

        Real lambda_in = (length(aux_in) / abs(aux_in.z) - 1) / 2;
        Real lambda_out = (length(aux_out) / abs(aux_out.z) - 1) / 2;

        masking = 1 / ((1 + lambda_in) * (1 + lambda_out));
        // masking = smith_masking_gtr2(dir_local_in, roughness) * smith_masking_gtr2(dir_local_out, roughness);
    }

    Real f = fresnel * distribution * masking / (4 * dot_ni);

    return make_const_spectrum(f);
}

Real pdf_sample_bsdf_op::operator()(const DisneyClearcoat &bsdf) const {
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) < 0) {
        frame = -frame;
    }

    Vector3 dir_micronormal = normalize(dir_in + dir_out);

    Real clearcoat_gloss = eval(bsdf.clearcoat_gloss, vertex.uv, vertex.uv_screen_size, texture_pool);
    
    Real dot_nmicron = abs(dot(frame.n, dir_micronormal));
    Real dot_microno = abs(dot(dir_micronormal, dir_out));

    Real distribution;
    {
        Real alpha = (1 - clearcoat_gloss) * Real(0.1) + clearcoat_gloss * Real(0.001);
        Real alpha2 = alpha * alpha;
        
        distribution = (alpha2 - 1) / (c_PI * log(alpha2) * (1 + (alpha2 - 1) * dot_nmicron * dot_nmicron));
    }

    Real pdf = distribution * dot_nmicron / (4 * dot_microno);

    return pdf;
}

std::optional<BSDFSampleRecord>
        sample_bsdf_op::operator()(const DisneyClearcoat &bsdf) const {
    if (dot(vertex.geometric_normal, dir_in) < 0) {
        // No light below the surface
        return std::nullopt;
    }
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) < 0) {
        frame = -frame;
    }

    Real clearcoat_gloss = eval(bsdf.clearcoat_gloss, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real alpha = (1 - clearcoat_gloss) * Real(0.1) + clearcoat_gloss * Real(0.001);

    Vector3 dir_local_micronormal;
    {
        Real cos_elevation = sqrt((1 - pow(alpha * alpha, 1 - rnd_param_uv.x)) / (1 - alpha * alpha));
        Real sin_elevation = sqrt(1 - cos_elevation * cos_elevation);

        Real cos_azimuth = cos(c_TWOPI * rnd_param_uv.y);
        // bahahahaha :joy: do not forget the sign of the squareroot, i just wasted 5 hours rederiving
        // all the importance sampling math and reading through what feels like heatz's entire collection
        // of works to figure out the bug in the code. only to realize that i messed up the basic 9th
        // grade trig identity sin^2 + cos^2 = 1, because taking squareroots creates ambiguity. :joy:
        Real sin_azimuth = sqrt(1 - cos_azimuth * cos_azimuth) * ((rnd_param_uv.y < Real(0.5)) * 2 - 1);

        dir_local_micronormal = { sin_elevation * cos_azimuth, sin_elevation * sin_azimuth, cos_elevation };
    }

    Vector3 dir_micronormal = to_world(frame, dir_local_micronormal);
    Vector3 dir_out = -dir_in + 2 * dot(dir_in, dir_micronormal) * dir_micronormal;

    if (dot(vertex.geometric_normal, dir_out) < 0) {
        return std::nullopt;
    }

    return BSDFSampleRecord { dir_out, 0, 1 };
}

TextureSpectrum get_texture_op::operator()(const DisneyClearcoat &) const {
    return make_constant_spectrum_texture(make_zero_spectrum());
}
