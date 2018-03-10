#pragma once

#include <afxwin.h>

#include <omp.h>

#include <list>

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
        double u1, u2, ua, q;
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
            1, 1, 1, 0.01
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

    struct particle
    {
        geom::point2d_t x, v;
        bool alive;
    };

    struct mesh_data
    {
        plot::world_t::ptr_t world;
        util::ptr_t < geom::mesh > mesh;
        util::ptr_t < std::vector < particle > > data;
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
        pd.world = plot::world_t::create();
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

            for each (auto & pt in *m.data)
            {
                point_painter.draw_at(dc, vp, pt.x);
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
            { -20 * p.w, -20 * p.h }, { 20 * p.w, -20 * p.h }, { 20 * p.w, 20 * p.h }, { -20 * p.w, 20 * p.h },
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
            md.mesh->add(geom::point2d_t(- 19 * p.w, j * p.dy), material::ext);
            md.mesh->add(geom::point2d_t(- 19 * p.w, j * p.dy), material::ext);
        }

        md.mesh->add(geom::point2d_t(- 19 * p.w, - 19 * p.h), material::ext);
        md.mesh->add(geom::point2d_t(- 19 * p.w, + 19 * p.h), material::ext);
        md.mesh->add(geom::point2d_t(19 * p.w, - 19 * p.h), material::ext);
        md.mesh->add(geom::point2d_t(19 * p.w, + 19 * p.h), material::ext);

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

        md.data->clear();
    }

    inline mesh_data make_system_data(const parameters & p)
    {
        mesh_data md;

        md.world = plot::world_t::create();

        md.data = util::create < std::vector < particle > > ();

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
            #pragma omp parallel for
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

    /* *************************************************
       ******** Finite Element Method (Galerkin) *******
       ************************************************* */

    /* adaptation of learning/ss-9-poisson-ag */
    class finel_galerkin
    {
    private:
        struct plane { double a, b, c; };
        struct sparse_matrix
        {
            std::vector < std::vector < std::pair < size_t, double > > > matrix;
        };
    private:
        util::ptr_t < geom::mesh > m;
        /* here and below: a x = c */
        /* sparse matrix dim(a) = dim(x) * dim(x) */
        sparse_matrix a;
        /* dim(c) = dim(x) */
        std::vector < double > c;
        std::vector < double > x;
        std::vector < double > x1;
        /* mapping: mesh vertice -> variable */
        std::vector < geom::mesh::idx_t > vars;
        const double accuracy_goal;
        const size_t iters;
        const parameters & p;
        const std::vector < double > & charges;
    public:
        finel_galerkin(const parameters & p,
                       util::ptr_t < geom::mesh > m,
                       const std::vector < double > & charges,
                       double accuracy_goal,
                       size_t iters = 100)
            : p(p)
            , m(m)
            , iters(iters)
            , accuracy_goal(accuracy_goal)
            , charges(charges)
        {
        }
    public:
        void next(std::vector < double > & r);
        void init();
        void update();
    private:
        void _next();
        bool _check_accuracy();
        plane _make_plane(geom::mesh::idx_t t,
                          geom::mesh::idx_t v) const;
        /* < grad phi_i | grad phi_k > */
        double _dot(geom::mesh::idx_t i,
                    geom::mesh::idx_t j) const;
        /* < rho_i | phi_k > */
        double _charge_dot(geom::mesh::idx_t i,
                           geom::mesh::idx_t j) const;
        bool _is_var(geom::mesh::idx_t v) const;
        double _area(geom::mesh::idx_t t) const;
        double _potential_of(geom::mesh::idx_t i) const;
    };

    inline bool finel_galerkin::_is_var(geom::mesh::idx_t v) const
    {
        return (m->flags_at(v) &
                (geom::mesh::superstruct | material::ext | material::grid)) == 0;
    }

    inline double finel_galerkin::_area(geom::mesh::idx_t t) const
    {
        auto p1 = m->point_at(m->triangles()[t].vertices[1]) -
            m->point_at(m->triangles()[t].vertices[0]);
        auto p2 = m->point_at(m->triangles()[t].vertices[2]) -
            m->point_at(m->triangles()[t].vertices[0]);
        return std::abs(p1.x * p2.y - p1.y * p2.x) / 2;
    }

    inline double finel_galerkin::_potential_of(geom::mesh::idx_t i) const
    {
        if (m->flags_at(i) & material::cathode)
            return p.ua;
        else if (m->flags_at(i) & material::grid1)
            return p.u1;
        else if (m->flags_at(i) & material::grid2)
            return p.u2;
        return 0;
    }

    inline finel_galerkin::plane finel_galerkin::_make_plane(
        geom::mesh::idx_t t, geom::mesh::idx_t v) const
    {
        auto & ti = m->triangles()[t];
        auto & p1 = m->point_at(ti.vertices[0]);
        auto & p2 = m->point_at(ti.vertices[1]);
        auto & p3 = m->point_at(ti.vertices[2]);
        double d = p1.x * (p2.y - p3.y) +
                   p2.x * (p3.y - p1.y) +
                   p3.x * (p1.y - p2.y);
        if (ti.vertices[0] == v)
        {
            return
            {
                (p2.y - p3.y) / d,
                (p3.x - p2.x) / d,
                (p2.x * p3.y - p3.x * p2.y) / d
            };
        }
        else if (ti.vertices[1] == v)
        {
            return
            {
                (p3.y - p1.y) / d,
                (p1.x - p3.x) / d,
                (p3.x * p1.y - p1.x * p3.y) / d
            };
        }
        else // if (ti.vertices[2] == v)
        {
            return
            {
                (p1.y - p2.y) / d,
                (p2.x - p1.x) / d,
                (p1.x * p2.y - p2.x * p1.y) / d
            };
        }
    }

    inline double finel_galerkin::_dot(
        geom::mesh::idx_t i, geom::mesh::idx_t j) const
    {
        double r = 0;
        auto & nt = m->vertices()[i].neighbor_triangles;
        for (auto it1 = nt.begin(); it1 != nt.end(); ++it1)
        {
            auto & ti = m->triangles()[*it1];
            if ((ti.vertices[0] != j) &&
                (ti.vertices[1] != j) &&
                (ti.vertices[2] != j))
                continue;
            auto p1 = _make_plane(*it1, i);
            for (auto it2 = nt.begin(); it2 != nt.end(); ++it2)
            {
                if (*it1 != *it2)
                    continue;
                auto p2 = _make_plane(*it1, j);
                r += - _area(*it1) * (
                    p1.a * p2.a +
                    p1.b * p2.b
                );
            }
        }
        return r;
    }

    inline double finel_galerkin::_charge_dot(
        geom::mesh::idx_t i, geom::mesh::idx_t j) const
    {
        double r = 0;
        auto & nt = m->vertices()[i].neighbor_triangles;
        for (auto it1 = nt.begin(); it1 != nt.end(); ++it1)
        {
            auto & ti = m->triangles()[*it1];
            if ((ti.vertices[0] != j) &&
                (ti.vertices[1] != j) &&
                (ti.vertices[2] != j))
                continue;
            auto p1 = _make_plane(*it1, i);
            p1.a *= charges[i]; p1.b *= charges[i]; p1.c *= charges[i];
            for (auto it2 = nt.begin(); it2 != nt.end(); ++it2)
            {
                if (*it1 != *it2)
                    continue;
                auto p2 = _make_plane(*it1, j);
                /* move to a coordinate system in which
                   the triangle is rectangular:

                   (x, y) -> q0 + u (r0 - q0) + v (s0 - q0)

                   integral over triangle f(x,y)dxdy will
                   transform into integral |J| f(u,v)dudv,
                   where u = (0..1), v = (0..1-u), |J| is
                   Jacobian of the coordinate transformation.
                 */
                plane uv1, uv2;
                auto & q0 = m->point_at(ti.vertices[0]);
                auto & r0 = m->point_at(ti.vertices[1]);
                auto & s0 = m->point_at(ti.vertices[2]);
                uv1.c = p1.c + p1.a * s0.x + p1.b * s0.y;
                uv1.a = p1.a * (q0.x - s0.x) + p1.b * (q0.y - s0.y);
                uv1.b = p1.a * (r0.x - s0.x) + p1.b * (r0.y - s0.y);
                uv2.c = p2.c + p2.a * s0.x + p2.b * s0.y;
                uv2.a = p2.a * (q0.x - s0.x) + p2.b * (q0.y - s0.y);
                uv2.b = p2.a * (r0.x - s0.x) + p2.b * (r0.y - s0.y);
                double d = r0.x * (s0.y - q0.y) +
                           q0.x * (r0.y - s0.y) +
                           s0.x * (q0.y - r0.y);
                r += std::abs(d) / 24 * (uv2.a * uv1.b + uv1.a * uv2.b +
                               2 * uv1.b * uv2.b + 2 * uv1.a * uv2.a +
                               4 * uv2.a * uv1.c + 4 * uv2.b * uv1.c +
                               4 * uv1.b * uv2.c + 4 * uv1.a * uv2.c +
                               12 * uv1.c * uv2.c);
            }
        }
        return r;
    }

    inline void finel_galerkin::init()
    {
        if (!vars.empty()) return;

        vars.resize(m->vertices().size());
        geom::mesh::idx_t var = 0;
        for (geom::mesh::idx_t i = 0; i < m->vertices().size(); ++i)
        {
            if (!_is_var(i)) vars[i] = SIZE_T_MAX;
            else             vars[i] = var++;
        }

        x.resize(var);
        x1.resize(var);

        c.resize(var);
        a.matrix.resize(var);

        #pragma omp parallel for
        for (int j = 0; j < (int)m->vertices().size(); ++j)
        {
            if (vars[j] == SIZE_T_MAX) continue;
            #pragma omp parallel for firstprivate(j)
            for (int i = 0; i < (int)m->vertices().size(); ++i)
            {
                if (vars[i] != SIZE_T_MAX)
                {
                    auto d = _dot(i, j);
                    if (d != 0)
                    {
                        #pragma omp critical
                        a.matrix[vars[j]].emplace_back(vars[i], d);
                    }
                }
                else
                {
                    c[vars[j]] += - _potential_of(i) * _dot(i, j);
                }
                c[vars[j]] += _charge_dot(i, j);
            }
        }
    }

    inline void finel_galerkin::update()
    {
        #pragma omp parallel for
        for (int j = 0; j < (int)m->vertices().size(); ++j)
        {
            if (vars[j] == SIZE_T_MAX) continue;
            c[vars[j]] = 0;
            #pragma omp parallel for firstprivate(j)
            for (int i = 0; i < (int)m->vertices().size(); ++i)
            {
                if (vars[i] == SIZE_T_MAX)
                    c[vars[j]] += - _potential_of(i) * _dot(i, j);
                c[vars[j]] += _charge_dot(i, j);
            }
        }
    }

    inline void finel_galerkin::_next()
    {
        /* one iteration of Kaczmarz method
           without any accuracy checks */
        double s1, s2, t;
        for (size_t i = 0; i < a.matrix.size(); ++i)
        {
            s1 = s2 = 0;
            for (size_t k = 0; k < a.matrix[i].size(); ++k)
            {
                size_t j = std::get < 0 > (a.matrix[i][k]);
                double v = std::get < 1 > (a.matrix[i][k]);
                s1 += v * x[j];
                s2 += v * v;
            }
            t = (c[i] - s1) / s2;
            for (size_t k = 0; k < a.matrix[i].size(); ++k)
            {
                x[std::get < 0 > (a.matrix[i][k])] +=
                    std::get < 1 > (a.matrix[i][k]) * t;
            }
        }
    }

    inline bool finel_galerkin::_check_accuracy()
    {
        double s1 = 0, s2 = 0;

        #pragma omp parallel for reduction(+:s1,s2)
        for (int i = 0; i < (int)x.size(); ++i)
        {
            s1 += (x[i] - x1[i]) * (x[i] - x1[i]);
            s2 += x[i] * x[i];
        }

        return s1 <= accuracy_goal * s2;
    }

    inline void finel_galerkin::next(std::vector < double > & r)
    {
        x[0] = x1[0] = 0.5f;
        #pragma omp parallel for
        for (int i = 1; i < (int)x.size(); ++i)
            x[i] = x1[i] = 0;

        for (;;)
        {
            for (size_t i = 0; i < iters; ++i) _next();
            if (_check_accuracy()) break;
            x1 = x;
        };

        r.resize(m->vertices().size());
        #pragma omp parallel for
        for (int k = 0; k < (int)m->vertices().size(); ++k)
        {
            if (vars[k] != SIZE_T_MAX)
                r[k] = x[vars[k]];
            else
                r[k] = _potential_of(k);
        }
    }

    /* *************************************************
       ************ Particle-Grid Simulation ***********
       ************************************************* */

    class particle_particle
    {
    private:
        std::vector < particle > & particles;
        util::ptr_t < geom::mesh > m;
        const parameters & p;
        std::vector < double > areas;
    public:
        std::vector < double > charges;
        std::vector < double > potential;
        std::vector < geom::point2d_t > field;
        std::list < size_t > dead;
    public:
        particle_particle(const parameters & p,
                          util::ptr_t < geom::mesh > m,
                          std::vector < particle > & particles)
            : p(p)
            , m(std::move(m))
            , particles(particles)
        {
        }
    public:
        void init();
        double collect_charges();
        void generate_particles();
        void adjust_particles();
    private:

        double _area(geom::mesh::idx_t t) const
        {
            auto p1 = m->point_at(m->triangles()[t].vertices[1]) -
                m->point_at(m->triangles()[t].vertices[0]);
            auto p2 = m->point_at(m->triangles()[t].vertices[2]) -
                m->point_at(m->triangles()[t].vertices[0]);
            return std::abs(p1.x * p2.y - p1.y * p2.x) / 2;
        }

        bool _grad(geom::mesh::idx_t t, geom::point2d_t & p) const;

        void _calc_field();

        void _adjust_particle(particle & pt, const geom::point2d_t & n) const;
    };

    inline void particle_particle::init()
    {
        field.resize(m->triangles().size());
        potential.resize(m->vertices().size());
        charges.resize(m->vertices().size());
        areas.resize(m->triangles().size());

        #pragma omp parallel for
        for (int i0 = 0; i0 < (int)m->triangles().size(); ++i0)
        {
            geom::mesh::idx_t i(i0);
            areas[i] = _area(i);
        }
    }

    inline double particle_particle::collect_charges()
    {
        size_t cur = 0;
        #pragma omp parallel for
        for (int i0 = 0; i0 < (int)charges.size(); ++i0)
        {
            geom::mesh::idx_t i(i0);
            charges[i] = 0;
        }
        std::vector < std::list < size_t > > local_dead;
        #pragma omp parallel
        {
            const auto cur_thread = omp_get_thread_num();
            const auto all_threads = omp_get_num_threads();
            #pragma omp single
            {
                local_dead.resize(all_threads);
            }
            #pragma omp for
            for (int i = 0; i < (int)particles.size(); ++i)
            {
                if (!particles[i].alive) continue;
                if (particles[i].x.x >= p.w) ++cur;
                if ((particles[i].x.x >= p.w) ||
                    (particles[i].x.x <= 0) ||
                    (particles[i].x.y >= p.h / 2) ||
                    (particles[i].x.y <= -p.h / 2))
                {
                    particles[i].alive = false;
                    local_dead[cur_thread].push_back(i);
                    continue;
                }
                auto dc = m->find_triangle(particles[i].x);
                if (dc != SIZE_T_MAX)
                {
                    #pragma omp atomic
                    charges[m->triangles()[dc].vertices[0]] += p.q / areas[dc] / 3;
                    #pragma omp atomic
                    charges[m->triangles()[dc].vertices[1]] += p.q / areas[dc] / 3;
                    #pragma omp atomic
                    charges[m->triangles()[dc].vertices[2]] += p.q / areas[dc] / 3;
                }
            }
        }
        for (size_t i = 0; i < local_dead.size(); ++i)
        {
            dead.splice(dead.end(), std::move(local_dead[i]));
        }
        return cur * p.q;
    }

    inline void particle_particle::generate_particles()
    {
        size_t n = (size_t)std::floor(p.i0 / p.dt);
        for (size_t i = 0; i < n; ++i)
        {
            model::particle p0
            {
                geom::point2d_t{ (rand() / (RAND_MAX + 1.) + 1) * p.w / 100,
                                 (rand() / (RAND_MAX + 1.) - 0.5) * p.h },
                geom::point2d_t{ (1 + rand() / (RAND_MAX + 1.) * 0.1) * p.v0, 0 },
                true
            };
            if (!dead.empty())
            {
                particles[dead.back()] = p0;
                dead.pop_back();
            }
            else
            {
                particles.push_back(p0);
            }
        }
    }

    inline bool particle_particle::_grad(
        geom::mesh::idx_t t, geom::point2d_t & p) const
    {
        auto & ti = m->triangles()[t];
        auto & p1 = m->point_at(ti.vertices[0]);
        auto & p2 = m->point_at(ti.vertices[1]);
        auto & p3 = m->point_at(ti.vertices[2]);

        double d0 = (p2.x - p3.x) * p1.y +
                    (p3.x - p1.x) * p2.y +
                    (p1.x - p2.x) * p3.y;

        double n1 = (p3.y - p2.y) * potential[ti.vertices[0]] +
                    (p1.y - p3.y) * potential[ti.vertices[1]] +
                    (p2.y - p1.y) * potential[ti.vertices[2]];
        double n2 = (p2.x - p3.x) * potential[ti.vertices[0]] +
                    (p3.x - p1.x) * potential[ti.vertices[1]] +
                    (p1.x - p2.x) * potential[ti.vertices[2]];

        p.x = n1 / d0;
        p.y = n2 / d0;

        if (!isfinite(p.x) || !isfinite(p.y))
            return false;

        return true;
    }

    inline void particle_particle::_calc_field()
    {
        geom::point2d_t n;
        #pragma omp parallel for private(n)
        for (int i0 = 0; i0 < (int) m->triangles().size(); ++i0)
        {
            geom::mesh::idx_t i(i0);
            if (_grad(i, n))
                field[i] = n;
            else
                field[i] = {};
        }
    }

    inline void particle_particle::adjust_particles()
    {
        _calc_field();
        #pragma omp parallel for
        for (int i = 0; i < (int)particles.size(); ++i)
        {
            if (!particles[i].alive) continue;
            auto t = m->find_triangle(particles[i].x);
            if (t == SIZE_T_MAX)
                continue;
            _adjust_particle(particles[i], field[t]);
        }
    }

    inline void particle_particle::_adjust_particle(
        particle & pt, const geom::point2d_t & n) const
    {
        pt.v += geom::point2d_t{ n.x * p.dt, n.y * p.dt };
        pt.x += geom::point2d_t{ pt.v.x * p.dt, pt.v.y * p.dt };
    }
}
