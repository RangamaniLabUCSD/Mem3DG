#include <iostream>
#include "ddgsolver/icosphere.h"
#include "geometrycentral/utilities/vector3.h"
#include <math.h>

using Vector3 = geometrycentral::Vector3;

size_t getMidPoint(size_t t1, size_t t2, std::vector<Vector3>& coords) {
	Vector3 p1 = coords[t1];
	Vector3 p2 = coords[t2];
	Vector3 pm = (p1 + p2) / 2;
	pm = pm.normalize();
	size_t i = coords.size();
	coords.push_back(pm);
	return i;
}

void icosphere(std::vector<Vector3>& coords, std::vector<std::vector<std::size_t>>& polygons, int n) {
	// initialize the vertex coordinate
	double t = (1 + sqrt(5)) / 2;

	coords.push_back(Vector3{ -1, t, 0 });
	coords.push_back(Vector3{ 1, t, 0 });
	coords.push_back(Vector3{ -1, -t, 0 });
	coords.push_back(Vector3{ 1, -t, 0 });
	coords.push_back(Vector3{ 0, -1, t });
	coords.push_back(Vector3{ 0, 1, t });
	coords.push_back(Vector3{ 0, -1, -t });
	coords.push_back(Vector3{ 0, 1, -t });
	coords.push_back(Vector3{ t, 0, -1 });
	coords.push_back(Vector3{ t, 0, 1 });
	coords.push_back(Vector3{ -t, 0, -1 });
	coords.push_back(Vector3{ -t, 0, 1 });
	
	
	for (int i = 0; i < coords.size(); ++i) {
		coords[i] = coords[i].normalize();
		//std::cout << "normalized" <<  coords[i] << std::endl;
	}
	/*
	for (Vector3 v : coords) {
		std::cout << v << std::endl;
		std::cout << "radius squared" << v[0] * v[0] + v[1] * v[1] + v[2] * v[2] << std::endl;
	}
	*/
	
	/*
	std::cout << "coords 2 1 is : " << coords[2][1] << std::endl;
	// if 2,4 it actually goes to 3,1
	std::cout << "coords 2 4 is : " << coords[2][4] << std::endl;
	*/

	// initialize the face
	polygons.push_back(std::vector<std::size_t>{-1+1,  -1+12, -1 + 6		});
	polygons.push_back(std::vector<std::size_t>{-1 + 1, -1 + 6, -1 + 2				});
	polygons.push_back(std::vector<std::size_t>{-1 + 1, -1 + 2, -1 + 8			});
	polygons.push_back(std::vector<std::size_t>{-1 + 1, -1 + 8, -1 + 11			});
	polygons.push_back(std::vector<std::size_t>{-1 + 1, -1 + 11, -1 + 12		});
	polygons.push_back(std::vector<std::size_t>{-1 + 2, -1 + 6, -1 + 10			});
	polygons.push_back(std::vector<std::size_t>{-1 + 6, -1 + 12, -1 + 5		});
	polygons.push_back(std::vector<std::size_t>{-1 + 12, -1 + 11, -1 + 3			});
	polygons.push_back(std::vector<std::size_t>{-1 + 11, -1 + 8, -1 + 7		});
	polygons.push_back(std::vector<std::size_t>{-1 + 8, -1 + 2, -1 + 9				});
	polygons.push_back(std::vector<std::size_t>{-1 + 4, -1 + 10, -1 + 5		});
	polygons.push_back(std::vector<std::size_t>{-1 + 4, -1 + 5, -1 + 3				});
	polygons.push_back(std::vector<std::size_t>{-1 + 4, -1 + 3, -1 + 7				});
	polygons.push_back(std::vector<std::size_t>{-1 + 4, -1 + 7, -1 + 9					});
	polygons.push_back(std::vector<std::size_t>{-1 + 4, -1 + 9, -1 + 10			});
	polygons.push_back(std::vector<std::size_t>{-1 + 5, -1 + 10, -1 + 6				});
	polygons.push_back(std::vector<std::size_t>{-1 + 3, -1 + 5, -1 + 12			});
	polygons.push_back(std::vector<std::size_t>{-1 + 7, -1 + 3, -1 + 11				});
	polygons.push_back(std::vector<std::size_t>{-1 + 9, -1 + 7, -1 + 8				});
	polygons.push_back(std::vector<std::size_t>{-1 + 10, -1 + 9, -1 + 2				});
	


	for (size_t iter = 0; iter < n; ++iter) {
		std::vector<std::vector<std::size_t>> polygons_new;
		for (size_t f = 0; f < polygons.size(); ++f) {
			auto triangle = polygons[f];
			size_t a = getMidPoint(triangle[0], triangle[1], coords);
			size_t b = getMidPoint(triangle[1], triangle[2], coords);
			size_t c = getMidPoint(triangle[2], triangle[0], coords);

			std::vector<std::vector<std::size_t>> new_face =
			{
				{triangle[0],a,c},
				{triangle[1],b,a},
				{triangle[2],c,b},
				{a,b,c}
			};
			

			for (size_t i = 4 * f; i < 4 * f + 4; ++i) {
				polygons_new.push_back(new_face[i - (4*f)]);
			}
			
		}
		polygons = polygons_new;
	}

	


}

