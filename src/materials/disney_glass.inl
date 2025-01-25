#include "../microfacet.h"

Spectrum eval_op::operator()(const DisneyGlass &bsdf) const {
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) * dot(vertex.geometric_normal, dir_in) < 0) {
        frame = -frame;
    }

    Spectrum base_color = eval(bsdf.base_color, vertex.uv, vertex.uv_screen_size, texture_pool);
    // Clamp roughness to avoid numerical issues.
    Real roughness = std::clamp(eval(bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool), Real(0.01), Real(1));
    Real anisotropic = eval(bsdf.anisotropic, vertex.uv, vertex.uv_screen_size, texture_pool);

    Real aspect = sqrt(1 - Real(0.9) * anisotropic);
    Real alpha_x = max(Real(0.0001), roughness * roughness / aspect);
    Real alpha_y = max(Real(0.0001), roughness * roughness * aspect);

    bool is_reflection = dot(vertex.geometric_normal, dir_in) * dot(vertex.geometric_normal, dir_out) > 0;

    Real ior_in = dot(vertex.geometric_normal, dir_in) < 0 ? bsdf.eta : 1;

    // Note ior_other is not necessarily the IOR of the light direction's medium. Only in
    // the case of refraction do they match. For computing the micronormal and BSDF, the
    // IOR of the light direction's medium is useful. But for computing the Fresnel
    // term, we always want the IOR of the mediums on *both* sides of the intersection.

    // The use of "other" instead of "transmitted" for naming is intentional. While the
    // Fresnel term is symmetric w.r.t. Snell's law, it feels better to call it the
    // "other" side. Even in the event of transmission, it is the view/in direction
    // that is the transmitted direction, physically speaking. The light's side is
    // not really the transmitted side.
    Real ior_other = bsdf.eta / ior_in;
    Real ior_out = is_reflection ? ior_in : ior_other;
    
    // Relative IOR in the literature refers to light:view, which is out:in from the
    // view's perspective (and lajolla uses the view's perspective for variable naming).
    Real ior_relative = ior_out / ior_in;

    // Note that in the case of reflection, the micronormal is the half vector (as expected)
    // It is oriented towards the lighter material (lower IOR), which should be air and
    // should generally correspond to the macronormal (except maybe in extreme cases? idk)
    Vector3 dir_micronormal = -normalize(ior_relative * dir_out + dir_in);

    Real dot_ni = abs(dot(frame.n, dir_in));
    Real dot_microni = abs(dot(dir_micronormal, dir_in));
    Real dot_microno = abs(dot(dir_micronormal, dir_out));

    if (length_squared(dir_micronormal) == 0) {
        return make_zero_spectrum();
    }

    Real fresnel_reflection;
    {
        Real ior_ratio = ior_in / ior_other;
        Real dot_micronother2 = 1 - ior_ratio * ior_ratio * (1 - dot_microni * dot_microni);
        if (dot_micronother2 < 0) {
            fresnel_reflection = 1;
        } else {
            Real dot_micronother = sqrt(dot_micronother2);

            Real R_s = (ior_in * dot_microni - ior_other * dot_micronother) / (ior_in * dot_microni + ior_other * dot_micronother);
            Real R_p = (ior_in * dot_micronother - ior_other * dot_microni) / (ior_in * dot_micronother + ior_other * dot_microni);
            fresnel_reflection = min(Real(1), ((R_s * R_s) + (R_p * R_p)) / 2);
        }
    }
    
    Real distribution;
    {
        Vector3 dir_local_half = to_local(frame, dir_micronormal);

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

    Spectrum f;
    {
        if (is_reflection) {
            f = base_color * fresnel_reflection * distribution * masking / (4 * dot_ni);
        } else {
            Real numer = (1 - fresnel_reflection) * distribution * masking * dot_microni * dot_microno;
            Real aux = dot(dir_micronormal, dir_in) + ior_relative * dot(dir_micronormal, dir_out);
            Real denom = dot_ni * aux * aux;

            // I'm extremely skeptical of needing a sqrt on the base color for refraction.
            // The argument that "it has two intersections, in + out" is not at all convincing.
            // I was also unable to find any references to this in the literature. But I could
            // be completely wrong here.
            f = sqrt(base_color) * numer / denom;
            // f = base_color * numer / denom;
        }
    }

    return f;
}

Real pdf_sample_bsdf_op::operator()(const DisneyGlass& bsdf) const
{
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) * dot(vertex.geometric_normal, dir_in) < 0) {
        frame = -frame;
    }

    // Clamp roughness to avoid numerical issues.
    Real roughness = std::clamp(eval(bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool), Real(0.01), Real(1));
    Real anisotropic = eval(bsdf.anisotropic, vertex.uv, vertex.uv_screen_size, texture_pool);

    Real aspect = sqrt(1 - Real(0.9) * anisotropic);
    Real alpha_x = max(Real(0.0001), roughness * roughness / aspect);
    Real alpha_y = max(Real(0.0001), roughness * roughness * aspect);

    bool is_reflection = dot(vertex.geometric_normal, dir_in) * dot(vertex.geometric_normal, dir_out) > 0;

    Real ior_in = dot(vertex.geometric_normal, dir_in) < 0 ? bsdf.eta : 1;
    Real ior_other = bsdf.eta / ior_in;
    Real ior_out = is_reflection ? ior_in : ior_other;
    Real ior_relative = ior_out / ior_in;

    Vector3 dir_micronormal = -normalize(ior_relative * dir_out + dir_in);

    Real dot_ni = abs(dot(frame.n, dir_in));
    Real dot_microni = abs(dot(dir_micronormal, dir_in));
    Real dot_microno = abs(dot(dir_micronormal, dir_out));

    Real fresnel_reflection;
    {
        Real ior_ratio = ior_in / ior_other;
        Real dot_micronother2 = 1 - ior_ratio * ior_ratio * (1 - dot_microni * dot_microni);
        if (dot_micronother2 < 0) {
            fresnel_reflection = 1;
        } else {
            Real dot_micronother = sqrt(dot_micronother2);

            Real R_s = (ior_in * dot_microni - ior_other * dot_micronother) / (ior_in * dot_microni + ior_other * dot_micronother);
            Real R_p = (ior_in * dot_micronother - ior_other * dot_microni) / (ior_in * dot_micronother + ior_other * dot_microni);
            fresnel_reflection = min(Real(1), ((R_s * R_s) + (R_p * R_p)) / 2);
        }
    }

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

    Real pdf;
    {
        if (is_reflection) {
            // Real visible_distribution = masking_in * dot_microni * distribution / dot_ni;
            // Real dh_dout = 1 / 4 * dot_microni;
            // Real pdf = fresnel * visible_distribution * dh_dout;

            // "Simplify" computation above; Can cancel out dot_microni to improve numerical stability
            pdf = fresnel_reflection * masking_in * distribution / (4 * dot_ni);
        } else {
            Real visible_distribution = masking_in * dot_microni * distribution / dot_ni;
            Real aux = dot(dir_micronormal, dir_in) + ior_relative * dot(dir_micronormal, dir_out);
            Real dh_dout = ior_relative * ior_relative * dot_microno / (aux * aux);
            pdf = (1 - fresnel_reflection) * visible_distribution * dh_dout;
        }
    }

    return pdf;
}

std::optional<BSDFSampleRecord> sample_bsdf_op::operator()(const DisneyGlass& bsdf) const
{
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) * dot(vertex.geometric_normal, dir_in) < 0) {
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
        Real sign = 2 * (dir_local_in.z >= 0) - 1;
        Frame frame_hemi_in { sign * normalize(Vector3 { alpha_x, alpha_y, Real(1) } * dir_local_in) };

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

        dir_local_micronormal = sign * normalize(vec_local_micronormal);
        // WARNING: sample_visible_normals does not support anisotropy
        // dir_local_micronormal = sample_visible_normals(dir_local_in, roughness * roughness, rnd_param_uv);
    }

    Vector3 dir_micronormal = to_world(frame, dir_local_micronormal);

    Real dot_microni = abs(dot(dir_micronormal, dir_in));

    Real ior_in = dot(vertex.geometric_normal, dir_in) < 0 ? bsdf.eta : 1;
    Real ior_other = bsdf.eta / ior_in;
    Real ior_ratio = ior_in / ior_other;

    Vector3 dir_other;
    Real fresnel_reflection;
    {
        Real dot_micronother2 = 1 - ior_ratio * ior_ratio * (1 - dot_microni * dot_microni);
        if (dot_micronother2 < 0) {
            fresnel_reflection = 1;
        } else {
            Real dot_micronother = sqrt(dot_micronother2);

            Real sign = 2 * (dot(dir_micronormal, dir_in) > 0) - 1;
            dir_other = sign * (ior_ratio * dot_microni - dot_micronother) * dir_micronormal - ior_ratio * dir_in;

            Real R_s = (ior_in * dot_microni - ior_other * dot_micronother) / (ior_in * dot_microni + ior_other * dot_micronother);
            Real R_p = (ior_in * dot_micronother - ior_other * dot_microni) / (ior_in * dot_micronother + ior_other * dot_microni);
            fresnel_reflection = min(Real(1), ((R_s * R_s) + (R_p * R_p)) / 2);
        }
    }

    Vector3 dir_out;
    Real ior_out;
    {
        if (rnd_param_w <= fresnel_reflection) {
            dir_out = -dir_in + 2 * dot(dir_in, dir_micronormal) * dir_micronormal;
            ior_out = ior_in;

            if (dot(vertex.geometric_normal, dir_out) * dot(vertex.geometric_normal, dir_in) < 0) {
                return std::nullopt;
            }
        } else {
            dir_out = dir_other;
            ior_out = ior_other;

            if (dot(vertex.geometric_normal, dir_out) * dot(vertex.geometric_normal, dir_in) > 0) {
                return std::nullopt;
            }
        }
    }

    return BSDFSampleRecord { dir_out, ior_out / ior_in, roughness };
}

TextureSpectrum get_texture_op::operator()(const DisneyGlass &bsdf) const {
    return bsdf.base_color;
}
