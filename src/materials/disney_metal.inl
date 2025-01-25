#include "../microfacet.h"

Spectrum eval_op::operator()(const DisneyMetal &bsdf) const {
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
    // Clamp roughness to avoid numerical issues.
    Real roughness = std::clamp(eval(bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool), Real(0.01), Real(1));
    Real anisotropic = eval(bsdf.anisotropic, vertex.uv, vertex.uv_screen_size, texture_pool);

    Vector3 dir_half = normalize(dir_in + dir_out);

    Real dot_ni = abs(dot(frame.n, dir_in));
    Real dot_ho = abs(dot(dir_half, dir_out));
    Real comp_ho = 1 - dot_ho;

    Real aspect = sqrt(1 - Real(0.9) * anisotropic);
    Real alpha_x = max(Real(0.0001), roughness * roughness / aspect);
    Real alpha_y = max(Real(0.0001), roughness * roughness * aspect);

    Spectrum fresnel = base_color + (make_const_spectrum(1) - base_color) * comp_ho * comp_ho * comp_ho * comp_ho * comp_ho;
    
    Real distribution;
    {
        Vector3 dir_local_half = to_local(frame, dir_half);

        Vector3 aux = dir_local_half / Vector3 { alpha_x, alpha_y, Real(1) };
        Real aux2 = dot(aux, aux);

        distribution = 1 / (c_PI * alpha_x * alpha_y * aux2 * aux2);
        // WARNING: GTR2 does not handle anisotropy
        // distribution = GTR2(dot(dir_half, frame.n), roughness);
    }

    Real masking;
    {
        Vector3 dir_local_in = to_local(frame, dir_in);
        Vector3 dir_local_out = to_local(frame, dir_out);

        Vector3 aux_in = Vector3 { alpha_x, alpha_y, Real(1) } * dir_local_in;
        Vector3 aux_out = Vector3 { alpha_x, alpha_y, Real(1) } * dir_local_out;

        Real lambda_in = (length(aux_in) / abs(aux_in.z) - 1) / 2;
        Real lambda_out = (length(aux_out) / abs(aux_out.z) - 1) / 2;

        masking = 1 / ((1 + lambda_in) * (1 + lambda_out));
        // WARNING: smith_masking_gtr2 does not handle anisotropy
        // masking = smith_masking_gtr2(dir_local_in, roughness) * smith_masking_gtr2(dir_local_out, roughness); 
    }

    Spectrum f = fresnel * distribution * masking / (4 * dot_ni);

    return f;
}

Real pdf_sample_bsdf_op::operator()(const DisneyMetal &bsdf) const {
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

    Vector3 dir_micronormal = normalize(dir_in + dir_out);

    // No light below the interpolated surface
    // No reversed microsurface
    // No light below the microsurface (not possible by sampling algorithm)
    if (dot(frame.n, dir_in) <= 0 || dot(frame.n, dir_micronormal) <= 0 /*|| dot(dir_micronormal, dir_in) <= 0*/) {
        return 0;
    }

    // Clamp roughness to avoid numerical issues.
    Real roughness = std::clamp(eval(bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool), Real(0.01), Real(1));
    Real anisotropic = eval(bsdf.anisotropic, vertex.uv, vertex.uv_screen_size, texture_pool);
    
    Real dot_ni = abs(dot(frame.n, dir_in));
    // Real dot_microni = abs(dot(dir_micronormal, dir_in));

    Real aspect = sqrt(1 - Real(0.9) * anisotropic);
    Real alpha_x = max(Real(0.0001), roughness * roughness / aspect);
    Real alpha_y = max(Real(0.0001), roughness * roughness * aspect);

    Real distribution;
    {
        Vector3 dir_local_micronormal = to_local(frame, dir_micronormal);

        Vector3 aux = dir_local_micronormal / Vector3 { alpha_x, alpha_y, Real(1) };
        Real aux2 = dot(aux, aux);

        distribution = 1 / (c_PI * alpha_x * alpha_y * aux2 * aux2);
        // WARNING: GTR2 does not handle anisotropy
        // distribution = GTR2(dot(dir_micronormal, frame.n), roughness);
    }

    Real masking_in;
    {
        Vector3 dir_local_in = to_local(frame, dir_in);

        Vector3 aux_in = Vector3 { alpha_x, alpha_y, Real(1) } * dir_local_in;
        Real lambda_in = (length(aux_in) / abs(aux_in.z) - 1) / 2;

        masking_in = 1 / (1 + lambda_in);
        // WARNING: smith_masking_gtr2 does not handle anisotropy
        // masking_in = smith_masking_gtr2(dir_local_in, roughness);
    }

    // Real visible_distribution = masking_in * dot_microni * distribution / dot_ni;
    // Real pdf = visible_distribution / (4 * dot_microni);

    // "Simplify" computation above; Can cancel out dot_microni to improve numerical stability
    Real pdf = masking_in * distribution / (4 * dot_ni);

    return pdf;
}

std::optional<BSDFSampleRecord>
        sample_bsdf_op::operator()(const DisneyMetal &bsdf) const {
    if (dot(vertex.geometric_normal, dir_in) < 0) {
        // No light below the surface
        return {};
    }
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) < 0) {
        frame = -frame;
    }

    // Clamp roughness to avoid numerical issues.
    Real roughness = std::clamp(eval(bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool), Real(0.01), Real(1));
    Real anisotropic = eval(bsdf.anisotropic, vertex.uv, vertex.uv_screen_size, texture_pool);

    Real aspect = sqrt(1 - Real(0.9) * anisotropic);
    Real alpha_x = max(Real(0.0001), roughness * roughness / aspect);
    Real alpha_y = max(Real(0.0001), roughness * roughness * aspect);

    Vector3 dir_local_in = to_local(frame, dir_in);

    Vector3 dir_local_micronormal;
    {

        Frame frame_hemi_in { normalize(Vector3 { alpha_x, alpha_y, Real(1) } * dir_local_in) };

        Real r = sqrt(rnd_param_uv.x);
        Real phi = c_TWOPI * rnd_param_uv.y;
        Real x = r * cos(phi);

        Real s = (1 + frame_hemi_in.n.z) / 2;
        Real y = (1 - s) * sqrt(1 - x * x) + s * r * sin(phi);

        Real z = sqrt(max(Real(0), 1 - x * x - y * y));

        Vector3 dir_in_micronormal = { x, y, z };
        Vector3 dir_hemi_micronormal = to_world(frame_hemi_in, dir_in_micronormal);

        Vector3 vec_local_micronormal = Vector3 { alpha_x, alpha_y, Real(1) } * dir_hemi_micronormal;
        vec_local_micronormal.z = max(Real(0), vec_local_micronormal.z);

        dir_local_micronormal = normalize(vec_local_micronormal);
        // WARNING: sample_visible_normals does not support anisotropy
        // dir_local_micronormal = sample_visible_normals(dir_local_in, roughness * roughness, rnd_param_uv);
    }

    Vector3 dir_micronormal = to_world(frame, dir_local_micronormal);
    Vector3 dir_out = -dir_in + 2 * dot(dir_in, dir_micronormal) * dir_micronormal;

    return BSDFSampleRecord { dir_out, 0, roughness };
}

TextureSpectrum get_texture_op::operator()(const DisneyMetal &bsdf) const {
    return bsdf.base_color;
}
