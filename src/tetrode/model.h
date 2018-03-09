#pragma once

#include <afxwin.h>

#include <util/common/geom/geom.h>
#include <util/common/plot/plot.h>
#include <util/common/plot/triangulation_drawable.h>
#include <util/common/plot/dirichlet_cell_drawable.h>

namespace model
{

    using points_t = std::vector < geom::point2d_t > ;

    struct parameters
    {
        // system params
        double w, h;
        double i0, v0;

        // geometry params
        double x1, x2;

        // other params
        double dx, dy, dt;
        size_t ndt;

        // material params
        double u1, u2, ua;
    };

    inline static parameters make_default_parameters()
    {
        return
        {
            // system params
            100, 10,
            10, 1,

            // geometry params
            35, 65,

            // other params
            10, 1, 0.01, 100,

            // material params
            1, 1, 1
        };
    }

    using material_t = geom::mesh::flags_t;

    namespace material
    {
        const material_t ext     = 0x1 << (0 + 10);
        const material_t cathode = 0x1 << (1 + 10);
        const material_t anode   = 0x1 << (2 + 10);
        const material_t grid1   = 0x1 << (3 + 10);
        const material_t grid2   = 0x1 << (4 + 10);
        const material_t grid    = cathode | anode | grid1 | grid2;
    };

    struct mesh_data
    {
        plot::world_t::ptr_t world;
        util::ptr_t < geom::mesh > mesh;
        util::ptr_t < points_t > data;
        plot::triangulation_drawable :: ptr_t triangulation_plot;
        plot::dirichlet_cell_drawable :: ptr_t dirichlet_cell_plot;
        plot::drawable::ptr_t system_plot;
        plot::drawable::ptr_t point_plot;
    };

    struct plot_data
    {
        plot::world_t::ptr_t world;
        util::ptr_t < points_t > data;
        plot::list_drawable < points_t > :: ptr_t plot;
    };

    struct multiplot_data
    {
        util::ptr_t < std::vector < points_t > > data;
        plot::multilist_drawable < points_t > :: ptr_t plot;
    };

    struct model_data
    {
        util::ptr_t < parameters > params;
        multiplot_data isoline_data;
        plot_data   func_data;
        mesh_data   system_data;
    };

    inline static plot_data make_plot_data
    (
        plot::palette::pen_ptr pen = plot::palette::pen(0xffffff)
    )
    {
        plot_data pd;
        pd.data = util::create < points_t > ();
        pd.plot = plot::list_drawable < points_t > :: create
        (
            plot::make_data_source(pd.data),
            nullptr, // no point painter
            pen
        );
        return pd;
    }

    inline static multiplot_data make_multiplot_data
    (
        plot::palette::pen_ptr pen = plot::palette::pen(0xffffff),
        plot::list_data_format data_format = plot::list_data_format::chain
    )
    {
        multiplot_data pd;
        pd.data = util::create < std::vector < points_t > > ();
        pd.plot = plot::multilist_drawable < points_t > :: create
        (
            plot::make_data_source(pd.data),
            nullptr, // no point painter
            pen
        );
        pd.plot->data_format = data_format;
        return pd;
    }

    inline static void adjust(const parameters & p,
                              plot::world_t & world)
    {
        world.xmin = 0;
        world.xmax = p.w;
        world.ymin = - (world.ymax = p.h / 2);
    }

    inline static plot::drawable::ptr_t make_root_drawable
    (
        plot::world_mapper_t world_mapper,
        std::vector < plot::drawable::ptr_t > layers
    )
    {
        using namespace plot;

        return viewporter::create(
            tick_drawable::create(
                layer_drawable::create(layers),
                const_n_tick_factory<axe::x>::create(
                    make_simple_tick_formatter(2, 5),
                    0,
                    5
                ),
                const_n_tick_factory<axe::y>::create(
                    make_simple_tick_formatter(2, 5),
                    0,
                    5
                ),
                palette::pen(RGB(80, 80, 80)),
                RGB(200, 200, 200)
            ),
            make_viewport_mapper(world_mapper)
        );
    }

    inline static plot::painter_t make_system_painter(const parameters & params,
                                                      mesh_data m)
    {
        return [&, m] (CDC & dc, const plot::viewport & vp)
        {
            auto point_brush = plot::palette::brush(0x3cff3c);

            auto point_painter = plot::custom_drawable(plot::circle_painter(3, point_brush));

            auto cathode_pen = plot::palette::pen(RGB(255, 255, 0), 8);
            auto anode_pen   = plot::palette::pen(RGB(255, 255, 0), 8);
            auto grid_pen    = plot::palette::pen(RGB(255, 255, 0), 4);

            dc.SelectObject(anode_pen.get());
            dc.MoveTo(vp.world_to_screen().xy({ 0, - params.h / 2 }));
            dc.LineTo(vp.world_to_screen().xy({ 0, + params.h / 2 }));
            dc.SelectObject(cathode_pen.get());
            dc.MoveTo(vp.world_to_screen().xy({ params.w, - params.h / 2 }));
            dc.LineTo(vp.world_to_screen().xy({ params.w, + params.h / 2 }));
            dc.SelectObject(grid_pen.get());
            dc.MoveTo(vp.world_to_screen().xy({ params.x1, - params.h / 2 }));
            dc.LineTo(vp.world_to_screen().xy({ params.x1, + params.h / 2 }));
            dc.MoveTo(vp.world_to_screen().xy({ params.x2, - params.h / 2 }));
            dc.LineTo(vp.world_to_screen().xy({ params.x2, + params.h / 2 }));

            for (size_t i = 0; i < m.data->size(); ++i)
            {
                point_painter.draw_at(dc, vp, m.data->at(i));
            }
        };
    }

    inline static plot::painter_t make_point_painter(const parameters & params,
                                                     mesh_data m)
    {
        return [&, m] (CDC & dc, const plot::viewport & vp)
        {
            auto metal_brush  = plot::palette::brush(0xacacac);

            auto point_painter = plot::custom_drawable(plot::circle_painter(2, metal_brush));

            for (geom::mesh::idx_t i = 0; i < m.mesh->vertices().size(); ++i)
            {
                point_painter.draw_at(dc, vp, m.mesh->point_at(i));
            }
        };
    }

    inline void update_system_data(const parameters & p, mesh_data & md)
    {
        std::vector < geom::point2d_t > super =
        {
            { -p.w, -20 * p.h }, { 2 * p.w, -20 * p.h }, { 2 * p.w, 20 * p.h }, { -p.w, 20 * p.h },
        };

        md.mesh->init(super);

        std::vector < geom::mesh::idx_t > vs;

        double grid_dy = p.dy / 2;
        size_t grid_n = size_t(std::floor(p.h / grid_dy + 1));

        for (size_t i = 0; i < grid_n; ++i)
        {
            md.mesh->add(geom::point2d_t(0,    - p.h / 2 + grid_dy * (double) i), material::anode);
            md.mesh->add(geom::point2d_t(p.w,  - p.h / 2 + grid_dy * (double) i), material::cathode);
            md.mesh->add(geom::point2d_t(p.x1, - p.h / 2 + grid_dy * (double) i), material::grid1);
            md.mesh->add(geom::point2d_t(p.x2, - p.h / 2 + grid_dy * (double) i), material::grid2);
        }

        size_t n = size_t(std::floor(p.w / p.dx + 1));
        size_t m = size_t(std::floor(p.h / p.dy + 1));
        for (size_t i = 1; i < n - 1; ++i)
        for (size_t j = 0; j < m; ++j)
        {
            md.mesh->add(geom::point2d_t(i * p.dx, -p.h / 2 + j * p.dy));
        }

        for (size_t i = 0; i < n; ++i)
        {
            md.mesh->add(geom::point2d_t(i * p.dx, -19 * p.h), material::ext);
            md.mesh->add(geom::point2d_t(i * p.dx, 19 * p.h), material::ext);
        }

        for (size_t j = 0; j < m; ++j)
        {
            md.mesh->add(geom::point2d_t(- p.w / 2, j * p.dy), material::ext);
            md.mesh->add(geom::point2d_t(- 3 * p.w / 2, j * p.dy), material::ext);
        }

        md.mesh->add(geom::point2d_t(- p.w / 2, - 19 * p.h), material::ext);
        md.mesh->add(geom::point2d_t(- p.w / 2, + 19 * p.h), material::ext);
        md.mesh->add(geom::point2d_t(3 * p.w / 2, - 19 * p.h), material::ext);
        md.mesh->add(geom::point2d_t(3 * p.w / 2, + 19 * p.h), material::ext);

        md.mesh->add(geom::point2d_t(0,    - 10 * p.h), material::anode);
        md.mesh->add(geom::point2d_t(p.w,  - 10 * p.h), material::cathode);
        md.mesh->add(geom::point2d_t(p.x1, - 10 * p.h), material::grid1);
        md.mesh->add(geom::point2d_t(p.x2, - 10 * p.h), material::grid2);
        md.mesh->add(geom::point2d_t(0,    + 10 * p.h), material::anode);
        md.mesh->add(geom::point2d_t(p.w,  + 10 * p.h), material::cathode);
        md.mesh->add(geom::point2d_t(p.x1, + 10 * p.h), material::grid1);
        md.mesh->add(geom::point2d_t(p.x2, + 10 * p.h), material::grid2);

        md.mesh->build_polygons = true;
        md.mesh->finish_mesh();
        md.mesh->build_neighborhood_tree();

        md.data = util::create < std::vector < geom::point2d_t > > ();
    }

    inline mesh_data make_system_data(const parameters & p)
    {
        mesh_data md;

        md.world = plot::world_t::create();

        md.mesh = util::create < geom::mesh > (false, false);

        update_system_data(p, md);

        md.triangulation_plot = plot::triangulation_drawable::create(
            plot::make_data_source(md.mesh), nullptr, plot::palette::pen(0xffffff));
        md.dirichlet_cell_plot = plot::dirichlet_cell_drawable::create(
            plot::make_data_source(md.mesh), nullptr, plot::palette::pen(RGB(255, 255, 0)));
        md.system_plot = plot::custom_drawable::create(
            make_system_painter(p, md));
        md.point_plot = plot::custom_drawable::create(
            make_point_painter(p, md));

        return md;
    }

    inline model_data make_model_data(const parameters & p = make_default_parameters())
    {
        model_data md;
        md.params = util::create < parameters > (p);
        md.func_data = make_plot_data(plot::palette::pen(0xffffff, 2));
        md.isoline_data = make_multiplot_data(plot::palette::pen(0xffffff, 2),
                                              plot::list_data_format::segment);
        md.system_data = make_system_data(*md.params);
        adjust(*md.params, *md.system_data.world);
        return md;
    }

    inline bool _get_val_intersection(double v1, double v2, double v, double & q)
    {
        double r = (v - v1) / (v2 - v1);
        if (isfinite(r) && (0 < r) && (r < 1))
        {
            q = r;
            return true;
        }
        return false;
    }

    inline void find_isolines(const geom::mesh & m,
                              const std::vector < double > & d,
                              double delta,
                              size_t max_lines,
                              std::vector < std::vector < geom::point2d_t > > & r)
    {
        r.clear();
        r.resize(max_lines * 2 + 1);
        for (geom::mesh::idx_t t = 0; t < m.triangles().size(); ++t)
        {
            auto & ti = m.triangles()[t];
            if (ti.flags & (geom::mesh::phantom | geom::mesh::superstruct)) continue;
            for (int l = - (int) max_lines; l <= (int) max_lines; ++l)
            {
                double val = l * delta;
                double q;
                size_t c = 0;
                if (_get_val_intersection(d[ti.vertices[0]], d[ti.vertices[1]], val, q))
                {
                    ++c;
                    r[(size_t) (l + max_lines)].push_back(
                        make_line_view(m.point_at(ti.vertices[0]),
                                       m.point_at(ti.vertices[1])).inner_point(q));
                }
                if (_get_val_intersection(d[ti.vertices[0]], d[ti.vertices[2]], val, q))
                {
                    ++c;
                    r[(size_t) (l + max_lines)].push_back(
                        make_line_view(m.point_at(ti.vertices[0]),
                                       m.point_at(ti.vertices[2])).inner_point(q));
                }
                if (_get_val_intersection(d[ti.vertices[1]], d[ti.vertices[2]], val, q))
                {
                    ++c;
                    r[(size_t) (l + max_lines)].push_back(
                        make_line_view(m.point_at(ti.vertices[1]),
                                       m.point_at(ti.vertices[2])).inner_point(q));
                }
                /* if occasionally added 1 or 3 points instead of 2;
                   may occur e.g. when the data contains !isfinite(data) */
                if ((c % 2) == 1)
                {
                    r[(size_t) (l + max_lines)].pop_back();
                }
            }
        }
    }
}
