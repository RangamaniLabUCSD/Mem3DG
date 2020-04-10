#include <iostream>
#include "ddgsolver/icosphere.h"
#include "geometrycentral/utilities/vector3.h"

using Vector3 = geometrycentral::Vector3;


void icosphere(std::vector<Vector3>& coords, std::vector<std::vector<std::size_t>>& polygons) {
	// initialize the vertex coordinate
	coords.push_back(Vector3{ 0.0f, 0.0f, 2.0f });
	coords.push_back(Vector3{ 1.788854f, 0.000000f, 0.894427f });
	coords.push_back(Vector3{ 0.552786f, 1.701302f, 0.894427f });
	coords.push_back(Vector3{ -1.447214f, 1.051462f, 0.894427f });
	coords.push_back(Vector3{ -1.447214f, -1.051462f, 0.894427f });
	coords.push_back(Vector3{ 0.552786f, -1.701302f, 0.894427f });
	coords.push_back(Vector3{ 1.447214f, 1.051462f, -0.894427f });
	coords.push_back(Vector3{ -0.552786f, 1.701302f, -0.894427f });
	coords.push_back(Vector3{ -1.788854f, 0.000000f, -0.894427f });
	coords.push_back(Vector3{ -0.552786f, -1.701302f, -0.894427f });
	coords.push_back(Vector3{ 1.447214f, -1.051462f, -0.894427f });
	coords.push_back(Vector3{ 0.0f , 0.0f,  -2.0f });
	for (Vector3 v : coords) {
		std::cout << v << std::endl;
		//std::cout << "radius squared" << v[0] * v[0] + v[1] * v[1] + v[2] * v[2] << std::endl;
	}
	/*
	std::cout << "coords 2 1 is : " << coords[2][1] << std::endl;
	// if 2,4 it actually goes to 3,1
	std::cout << "coords 2 4 is : " << coords[2][4] << std::endl;
	*/

	// initialize the face
	polygons.push_back(std::vector<std::size_t>{2, 0, 1});
	polygons.push_back(std::vector<std::size_t>{3, 0, 2	});
	polygons.push_back(std::vector<std::size_t>{4, 0, 3});
	polygons.push_back(std::vector<std::size_t>{5, 0, 4	});

	polygons.push_back(std::vector<std::size_t>{1, 0, 5});
	polygons.push_back(std::vector<std::size_t>{2, 1, 6	});
	polygons.push_back(std::vector<std::size_t>{7, 2, 6});
	polygons.push_back(std::vector<std::size_t>{3, 2, 7	});

	polygons.push_back(std::vector<std::size_t>{8, 3, 7});
	polygons.push_back(std::vector<std::size_t>{4, 3, 8	});
	polygons.push_back(std::vector<std::size_t>{9, 4, 8});
	polygons.push_back(std::vector<std::size_t>{5, 4, 9	});

	polygons.push_back(std::vector<std::size_t>{10, 5, 9});
	polygons.push_back(std::vector<std::size_t>{6, 1, 10	});
	polygons.push_back(std::vector<std::size_t>{1, 5, 10});
	polygons.push_back(std::vector<std::size_t>{6, 11, 7	});

	polygons.push_back(std::vector<std::size_t>{7, 11, 8});
	polygons.push_back(std::vector<std::size_t>{8, 11, 9	});
	polygons.push_back(std::vector<std::size_t>{9, 11, 10});
	polygons.push_back(std::vector<std::size_t>{10, 11, 6	});

}