#pragma once
#include <iostream>
#include <OpenMesh\Core\Mesh\PolyConnectivity.hh>
#include <OpenMesh\Core\Mesh\TriConnectivity.hh>
//��֪�������������Ǹ�ɶ�õ�?
//������ν���һ����ʼ��������Ҫ������Щͷ�ļ��أ�
#include <OpenMesh\Core\IO\MeshIO.hh>
#include <OpenMesh\Core\Mesh\TriMesh_ArrayKernelT.hh>
//#include <Include\OpenMesh\Core\Mesh\PolyMesh_ArrayKernelT.hh>//added by trytotry
//#include "MeshViewer\MeshDefinition.h"

typedef OpenMesh::TriMesh_ArrayKernelT<>  Mesh;
//typedef OpenMesh::PolyMesh_ArrayKernelT<> MyMesh;

class Remesh
{
public:
	Remesh();
	Remesh(Mesh *my_mesh_);
	~Remesh();

public:
	void initial_Para();

public:
	void split_long_edges();
	void collapse_short_edges();
	bool is_valences_ok(const OpenMesh::EdgeHandle &);
	void equalize_valences();
	void tangential_relaxation();

public:
	void remeshing();

private:
	Mesh				*my_mesh;
	double              ave_length_of_e;
	double              high_e;
	double              low_e;
};

