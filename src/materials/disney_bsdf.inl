#include "../microfacet.h"

Spectrum eval_op::operator()(const DisneyBSDF &bsdf) const {
    const auto modified_metal_bsdf_eval = [&]() -> Spectrum {
        // Flip the shading frame if it is inconsistent with the geometry normal
        Frame frame = vertex.shading_frame;
        if (dot(frame.n, dir_in) < 0) {
            frame = -frame;
        }

        Spectrum base_color = eval(bsdf.base_color, vertex.uv, vertex.uv_screen_size, texture_pool);
        // Clamp roughness to avoid numerical issues.
        Real roughness = std::clamp(eval(bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool), Real(0.01), Real(1));
        Real anisotropic = eval(bsdf.anisotropic, vertex.uv, vertex.uv_screen_size, texture_pool);

        Real specular = eval(bsdf.specular, vertex.uv, vertex.uv_screen_size, texture_pool);
        Real specular_tint = eval(bsdf.specular_tint, vertex.uv, vertex.uv_screen_size, texture_pool);
        Real metallic = eval(bsdf.metallic, vertex.uv, vertex.uv_screen_size, texture_pool);

        Vector3 dir_half = normalize(dir_in + dir_out);

        Real dot_ni = abs(dot(frame.n, dir_in));
        Real dot_ho = abs(dot(dir_half, dir_out));
        Real comp_ho = 1 - dot_ho;

        Real aspect = sqrt(1 - Real(0.9) * anisotropic);
        Real alpha_x = max(Real(0.0001), roughness * roughness / aspect);
        Real alpha_y = max(Real(0.0001), roughness * roughness * aspect);

        Spectrum fresnel;
        {
            Real l = luminance(base_color);
            Spectrum C_tint = l > 0 ? base_color / l : make_const_spectrum(1);
            Spectrum K_s = make_const_spectrum(1 - specular_tint) + specular_tint * C_tint;
            Real R_0 = (bsdf.eta - 1) * (bsdf.eta - 1) / ((bsdf.eta + 1) * (bsdf.eta + 1));
            Spectrum C_0 = specular * R_0 * (1 - metallic) * K_s + metallic * base_color;
            fresnel = C_0 + (make_const_spectrum(1) - C_0) * comp_ho * comp_ho * comp_ho * comp_ho * comp_ho;
        }

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
    };

    Real clearcoat = eval(bsdf.clearcoat, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real metallic = eval(bsdf.metallic, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real sheen = eval(bsdf.sheen, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real specular_transmission = eval(bsdf.specular_transmission, vertex.uv, vertex.uv_screen_size, texture_pool);

    Real diffuse_weight = (1 - specular_transmission) * (1 - metallic);
    Real sheen_weight = (1 - metallic) * sheen;
    // THIS IS INCORRECT. IT IS FOR THE ASSIGNMENT. The (modified metal weight)
    // is already computed by the modified fresnel, which lerps between two
    // (dielectric specular) and (metal) bsdfs. Having it here too is double
    // counting.
    Real metal_weight = metallic + (1 - specular_transmission) * (1 - metallic);
    // Real metal_weight = 1;
    Real clearcoat_weight = clearcoat / 4;
    Real glass_weight = (1 - metallic) * specular_transmission;

    if (dot(vertex.geometric_normal, dir_in) <= 0) {
        return glass_weight * (*this)(DisneyGlass { bsdf.base_color, bsdf.roughness, bsdf.anisotropic, bsdf.eta });
    }

    Spectrum f_diffuse = diffuse_weight * (*this)(DisneyDiffuse { bsdf.base_color, bsdf.roughness, bsdf.subsurface });
    Spectrum f_sheen = sheen_weight * (*this)(DisneySheen { bsdf.base_color, bsdf.sheen_tint });
    Spectrum f_metal = metal_weight * modified_metal_bsdf_eval();
    Spectrum f_clearcoat = clearcoat_weight * (*this)(DisneyClearcoat { bsdf.clearcoat_gloss });
    Spectrum f_glass = glass_weight * (*this)(DisneyGlass { bsdf.base_color, bsdf.roughness, bsdf.anisotropic, bsdf.eta });

    return f_diffuse + f_sheen + f_metal + f_clearcoat + f_glass;
}

Real pdf_sample_bsdf_op::operator()(const DisneyBSDF& bsdf) const {
    if (dot(vertex.geometric_normal, dir_in) <= 0) {
        return (*this)(DisneyGlass { bsdf.base_color, bsdf.roughness, bsdf.anisotropic, bsdf.eta });
    }

    Real clearcoat = eval(bsdf.clearcoat, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real metallic = eval(bsdf.metallic, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real specular_transmission = eval(bsdf.specular_transmission, vertex.uv, vertex.uv_screen_size, texture_pool);

    Real diffuse_weight = (1 - specular_transmission) * (1 - metallic);
    // The modified metal BSDF is actually a combined (metallic) + (specular dielectric),
    // so we have the (metallic weight) + (specular dielectric weight). The
    // (specular dielectric weight) is the same as the (diffuse weight).
    Real metal_weight = metallic + (1 - specular_transmission) * (1 - metallic);
    // Due to the incorrect double counting of the combined metal weight (see eval() above),
    // we square this for importance sampling.
    metal_weight *= metal_weight;
    Real clearcoat_weight = clearcoat / 4;
    Real glass_weight = (1 - metallic) * specular_transmission;
    Real total_weight = diffuse_weight + metal_weight + clearcoat_weight + glass_weight;

    if (total_weight == 0) {
        return (*this)(DisneyMetal { bsdf.base_color, bsdf.roughness, bsdf.anisotropic });
    }

    Real pr_diffuse = diffuse_weight / total_weight;
    Real pr_metal = metal_weight / total_weight;
    Real pr_clearcoat = clearcoat_weight / total_weight;
    Real pr_glass = glass_weight / total_weight;

    Real pdf_diffuse = pr_diffuse * (*this)(DisneyDiffuse { bsdf.base_color, bsdf.roughness, bsdf.subsurface });
    Real pdf_metal = pr_metal * (*this)(DisneyMetal { bsdf.base_color, bsdf.roughness, bsdf.anisotropic });
    Real pdf_clearcoat = pr_clearcoat * (*this)(DisneyClearcoat { bsdf.clearcoat_gloss });
    Real pdf_glass = pr_glass * (*this)(DisneyGlass { bsdf.base_color, bsdf.roughness, bsdf.anisotropic, bsdf.eta });

    Real pdf = max(Real(0), pdf_diffuse) + max(Real(0), pdf_metal) + max(Real(0), pdf_clearcoat) + max(Real(0), pdf_glass);

    return pdf;
}

std::optional<BSDFSampleRecord> sample_bsdf_op::operator()(const DisneyBSDF &bsdf) const {
    if (dot(vertex.geometric_normal, dir_in) <= 0) {
        return (*this)(DisneyGlass { bsdf.base_color, bsdf.roughness, bsdf.anisotropic, bsdf.eta });
    }

    Real clearcoat = eval(bsdf.clearcoat, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real metallic = eval(bsdf.metallic, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real specular_transmission = eval(bsdf.specular_transmission, vertex.uv, vertex.uv_screen_size, texture_pool);

    Real diffuse_weight = (1 - specular_transmission) * (1 - metallic);
    // The modified metal BSDF is actually a combined (metallic) + (specular dielectric),
    // so we have the (metallic weight) + (specular dielectric weight). The
    // (specular dielectric weight) is the same as the (diffuse weight).
    Real metal_weight = metallic + (1 - specular_transmission) * (1 - metallic);
    // Due to the incorrect double counting of the combined metal weight (see eval() above),
    // we square this for importance sampling.
    metal_weight *= metal_weight;
    Real clearcoat_weight = clearcoat / 4;
    Real glass_weight = (1 - metallic) * specular_transmission;
    Real total_weight = diffuse_weight + metal_weight + clearcoat_weight + glass_weight;

    if (total_weight == 0) {
        return (*this)(DisneyMetal { bsdf.base_color, bsdf.roughness, bsdf.anisotropic });
    }

    Real pr_diffuse = diffuse_weight / total_weight;
    Real pr_metal = metal_weight / total_weight;
    Real pr_clearcoat = clearcoat_weight / total_weight;
    Real pr_glass = glass_weight / total_weight;

    Real cdf_diffuse = 0 + pr_diffuse;
    Real cdf_metal = cdf_diffuse + pr_metal;
    Real cdf_clearcoat = cdf_metal + pr_clearcoat;
    Real cdf_glass = cdf_clearcoat + pr_glass;

    constexpr uint32_t offset = 1 << 24;
    double w;
    double w_glass = std::modf(rnd_param_w * offset, &w);
    w /= offset;

    if (w <= cdf_diffuse) {
        return (*this)(DisneyDiffuse { bsdf.base_color, bsdf.roughness, bsdf.subsurface });
    } else if (w <= cdf_metal) {
        return (*this)(DisneyMetal { bsdf.base_color, bsdf.roughness, bsdf.anisotropic });
    } else if (w <= cdf_clearcoat) {
        return (*this)(DisneyClearcoat { bsdf.clearcoat_gloss });
    } else if (w <= cdf_glass) {
        DisneyGlass bsdf_glass { bsdf.base_color, bsdf.roughness, bsdf.anisotropic, bsdf.eta };
        return sample_bsdf_op { dir_in, vertex, texture_pool, rnd_param_uv, w_glass, dir }(bsdf_glass);
    }

    return std::nullopt;
}

TextureSpectrum get_texture_op::operator()(const DisneyBSDF &bsdf) const {
    return bsdf.base_color;
}
