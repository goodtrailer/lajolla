Spectrum eval_op::operator()(const DisneyDiffuse &bsdf) const {
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
    Real roughness = eval(bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real subsurface = eval(bsdf.subsurface, vertex.uv, vertex.uv_screen_size, texture_pool);

    Vector3 dir_half = normalize(((dir_in + dir_out) / Real(2)));

    Real dot_ni = abs(dot(frame.n, dir_in));
    Real dot_no = abs(dot(frame.n, dir_out));
    Real dot_ho = dot(dir_half, dir_out);

    Real comp_ni = 1 - dot_ni;
    Real comp_no = 1 - dot_no;

    Spectrum f_base_diffuse;
    {
        Real fresnel_90 = Real(0.5) + 2 * roughness * dot_ho * dot_ho;
        Real fresnel_in = 1 + (fresnel_90 - 1) * comp_ni * comp_ni * comp_ni * comp_ni * comp_ni;
        Real fresnel_out = 1 + (fresnel_90 - 1) * comp_no * comp_no * comp_no * comp_no * comp_no;

        f_base_diffuse = base_color / c_PI * fresnel_in * fresnel_out * dot_no;
    }

    Spectrum f_subsurface;
    {
        Real fresnel_90 = roughness * dot_ho * dot_ho;
        Real fresnel_in = 1 + (fresnel_90 - 1) * comp_ni * comp_ni * comp_ni * comp_ni * comp_ni;
        Real fresnel_out = 1 + (fresnel_90 - 1) * comp_no * comp_no * comp_no * comp_no * comp_no;

        f_subsurface = Real(1.25) * base_color / c_PI * (fresnel_in * fresnel_out * (1 / (dot_ni + dot_no) - Real(0.5)) + Real(0.5)) * dot_no;
    }

    return (1 - subsurface) * f_base_diffuse + subsurface * f_subsurface;
}

Real pdf_sample_bsdf_op::operator()(const DisneyDiffuse &) const {
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

std::optional<BSDFSampleRecord> sample_bsdf_op::operator()(const DisneyDiffuse &) const {
    if (dot(vertex.geometric_normal, dir_in) < 0) {
        // No light below the surface
        return {};
    }
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) < 0) {
        frame = -frame;
    }

    return BSDFSampleRecord {
        to_world(frame, sample_cos_hemisphere(rnd_param_uv)),
        Real(0) /* eta */, Real(1) /* roughness */
    };
}

TextureSpectrum get_texture_op::operator()(const DisneyDiffuse &bsdf) const {
    return bsdf.base_color;
}
