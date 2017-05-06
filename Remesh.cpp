#include "Remesh.h"
#include <time.h>

Remesh::Remesh()
{
}


Remesh::~Remesh()
{
}

Remesh::Remesh(Mesh *my_mesh_)
{
	my_mesh = my_mesh_;
}

void Remesh::initial_Para()
{
	//compute the avearage length of edge of the mesh
	double sum = 0.;
	for (auto it = my_mesh->edges_begin(); it != my_mesh->edges_end(); ++it)
	{
		sum += my_mesh->calc_edge_length(it.handle());
	}
	ave_length_of_e = sum / my_mesh->n_edges();
	ave_length_of_e *= 0.95;

	/*std::cout << "ave_length_of_e : " << my_mesh->n_edges() << std::endl;
	std::cout << "ave_length_of_e : "<< ave_length_of_e << std::endl;*/

	high_e = 4 / 3 * ave_length_of_e;
	low_e = 4 / 5 * ave_length_of_e;
	//std::cout << "the number of edges : " << my_mesh->n_edges() << std::endl;
}

void Remesh::split_long_edges()
{
	int num = 0;
	for (auto it = my_mesh->edges_begin(); it != my_mesh->edges_end(); ++it)
		//it's some differences between "edges_begin()" and "edges_sbegin()"
		//what's the differences ?
		//it's about the mechanism of Openmesh; 
	{
		//std::cout << "The " << ++num << " times " << std::endl;
		if (my_mesh->calc_edge_length(it.handle()) > high_e)
		{
			OpenMesh::HalfedgeHandle _heh = my_mesh->halfedge_handle(it.handle(), 0);
			OpenMesh::VertexHandle point1 = my_mesh->from_vertex_handle(_heh);
			OpenMesh::VertexHandle point2 = my_mesh->to_vertex_handle(_heh);
			auto new_vert = (my_mesh->point(point1) + my_mesh->point(point2)) / 2;
			auto _vh = my_mesh->add_vertex(new_vert);
			//my_mesh->new_vertex()
			my_mesh->split(it.handle(),_vh);
		}
	}
}

void Remesh::collapse_short_edges()
{
	//细心 
	//仔细检查每个判断语句
	int num = 0;
	int num_of_collapse = 0;
	int num_of_boundary = 0;
	my_mesh->request_edge_status();
	my_mesh->request_face_status();
	my_mesh->request_vertex_status();
	for (auto it = my_mesh->edges_sbegin(); it != my_mesh->edges_end(); ++it)
	{
		//std::cout << "The " << ++num << " times " << std::endl;
		auto _heh = my_mesh->halfedge_handle(it.handle(), 0);
		if (my_mesh->is_boundary(it.handle()))
		{
			num_of_boundary++;
			continue;
		}
		
		//std::cout << "the error : " << my_mesh->calc_edge_length(it.handle()) - low_e<< std::endl;
		if (my_mesh->is_collapse_ok(_heh)  && my_mesh->calc_edge_length(it.handle()) < low_e)//这条语句需不需要加括号呢？
		{
			int flag = 1;
			//std::cout << "flag: " << flag << std::endl;
			OpenMesh::VertexHandle vert1 = my_mesh->from_vertex_handle(_heh);
			OpenMesh::VertexHandle vert2 = my_mesh->to_vertex_handle(_heh);
			auto point1 = my_mesh->point(vert1);
			auto point2 = my_mesh->point(vert2);

			auto new_vert = (point1 + point2) / 2;

			//第一个点的邻域
			for (auto it1 = my_mesh->vv_begin(vert1); it1 != my_mesh->vv_end(vert1); ++it1)
			{
				auto temp_point = my_mesh->point(it1);
				auto length = new_vert - temp_point;
				if (length.length() > high_e)
				{
					flag = 0;
					break;
				}
			}

			//第二个点的邻域
			for (auto it1 = my_mesh->vv_begin(vert2); it1 != my_mesh->vv_end(vert2); ++it1)
			{
				auto temp_point = my_mesh->point(it1);
				auto length = new_vert - temp_point;
				if (length.length() > high_e)
				{
					flag = 0;
					break;
				}
			}

			if (flag == 1)
			{
				//std::cout << "the number " << ++num_of_collapse << std::endl;
				my_mesh->set_point(vert2, new_vert);
				my_mesh->collapse(_heh);
			}

		}
	}
	/*std::cout << "num_of_boundary: " << num_of_boundary << std::endl;
	std::cout << "num_of_collapse: " << num_of_collapse << std::endl;*/
}

bool Remesh::is_valences_ok(const OpenMesh::EdgeHandle &_eh)
{
	auto _heh0 = my_mesh->halfedge_handle(_eh, 0); auto vert1 = my_mesh->to_vertex_handle(_heh0);
	auto _heh1 = my_mesh->halfedge_handle(_eh, 1); auto vert2 = my_mesh->to_vertex_handle(_heh1);
	auto _heh0_next = my_mesh->next_halfedge_handle(_heh0); auto vert3 = my_mesh->to_vertex_handle(_heh0_next);
	auto _heh1_next = my_mesh->next_halfedge_handle(_heh1); auto vert4 = my_mesh->to_vertex_handle(_heh1_next);

	bool b1 = my_mesh->is_boundary(vert1);
	bool b2 = my_mesh->is_boundary(vert2);
	bool b3 = my_mesh->is_boundary(vert3);
	bool b4 = my_mesh->is_boundary(vert4);

	int v1 = my_mesh->valence(vert1);
	int v2 = my_mesh->valence(vert2);
	int v3 = my_mesh->valence(vert3);
	int v4 = my_mesh->valence(vert4);

	int valence_ori = 0, valence_aft = 0;

	if (b1)
	{
		valence_ori += (v1 - 4)*(v1 - 4);
		valence_aft += (v1 - 5)*(v1 - 5);
	}
	else
	{
		valence_ori += (v1 - 6)*(v1 - 6);
		valence_aft += (v1 - 7)*(v3 - 7);
	}

	if (b2)
	{
		valence_ori += (v2 - 4)*(v2 - 4);
		valence_aft += (v2 - 5)*(v2 - 5);
	}
	else
	{
		valence_ori += (v2 - 6)*(v2 - 6);
		valence_aft += (v2 - 7)*(v2 - 7);
	}

	if (b3)
	{
		valence_ori += (v3 - 4)*(v3 - 4);
		valence_aft += (v3 - 3)*(v3 - 3);
	}
	else
	{
		valence_ori += (v3 - 6)*(v3 - 6);
		valence_aft += (v3 - 5)*(v3 - 5);
	}

	if (b4)
	{
		valence_ori += (v4 - 4)*(v4 - 4);
		valence_aft += (v4 - 3)*(v4 - 3);
	}
	else
	{
		valence_ori += (v4 - 6)*(v4 - 6);
		valence_aft += (v4 - 5)*(v4 - 5);
	}

	//仔细查看上面的代码有没有错误，有可能出错

	if (valence_ori > valence_aft)
		return true;

	return false;
}

void Remesh::equalize_valences()
{
	int num = 0;
	for (auto it = my_mesh->edges_sbegin(); it != my_mesh->edges_end(); ++it)
	{
		//std::cout << "The " << ++num << " times " << std::endl;
		if (!my_mesh->is_boundary(it.handle()) && my_mesh->is_flip_ok(it.handle()) && is_valences_ok(it.handle()))//看一下这些结构有没有错
			//flip_openmesh(it.handle(), *my_mesh);
			my_mesh->flip(it.handle());//is there any problem?
	}
}

void Remesh::tangential_relaxation()
{
	for (auto it = my_mesh->vertices_begin(); it != my_mesh->vertices_end(); ++it)
	{
		//compute the average vert
		OpenMesh::VectorT<float, 3> q; q[0] = 0.; q[1] = 0.; q[2] = 0.;
		for (auto it1 = my_mesh->vv_begin(*it); it1 != my_mesh->vv_end(*it); ++it1)
		{
			//auto p = my_mesh->point(*it1);
			q += my_mesh->point(*it1);
		}
		q /= my_mesh->valence(*it);

		my_mesh->request_vertex_normals();
		my_mesh->update_normals();

		auto normal = my_mesh->normal(*it);
		auto p = my_mesh->point(*it);
		double length = normal | (p - q);
		q = q + normal * length;
		my_mesh->set_point(*it, q);
	}
}

void Remesh::remeshing()
{
	//std::cout << "test remesh" << std::endl;
	time_t t_start, t_end;
	t_start = time(NULL);
	//for (size_t i = 0; i < 5; i++)
	//{
		initial_Para();
		split_long_edges();
		collapse_short_edges();
		equalize_valences();
		tangential_relaxation();
	//}
	t_end = time(NULL);

	printf("time: %.0f s\n", difftime(t_end, t_start));
}