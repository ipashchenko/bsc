#define CATCH_CONFIG_MAIN

#include "catch.hpp"
//
//
//TEST_CASE("Intersection of Cone with Ray", "[Cone]") {
//	// Just create Cone instance for subsequent usage
//	Vector3d origin(0.0, 0.0, 0.0);
//	Vector3d direction(0.0, 0.0, 1.0);
//	const double angle = pi/4.0;
//	const double scale = 10.0;
//	const double smallest_scale = 1E-05;
//	Cone cone(origin, direction, angle, scale);
//
//	SECTION("Two intersections") {
//		Vector3d ray_direction(0., 1., 0.);
//		Vector3d ray_origin(0., 1., 1.);
//		Vector3d point1(0.0, -1.0, 1.0);
//		Vector3d point2(0.0, 1.0, 1.0);
//		Ray ray(ray_origin, ray_direction);
//    std::list<Intersection> list_intersect = cone.hit(ray);
//		Intersection intersect = list_intersect.front();
//		std::pair<Vector3d,Vector3d> borders = intersect.get_path();
//		REQUIRE(list_intersect.size() == 1);
//		REQUIRE((borders.first - point1).isMuchSmallerThan(1E-06, 1E-06));
//		REQUIRE((borders.second - point2).isMuchSmallerThan(1E-06, 1E-06));
//	}
//
//	SECTION("No intersection") {
//		std::list<std::pair<Vector3d,Vector3d>> ray_directions;
//		double theta;
//		double r = 2.0;
//		double rr;
//		for (int j = 0; j < 8; ++j) {
//			theta = j*pi/8.;
//			for (int k = 0; k < 3; ++k) {
//				rr = k-1.0;
//				ray_directions.emplace_back(std::make_pair(Vector3d(r*cos(theta),
//				                                                    r*sin(theta),
//				                                                    rr),
//				                                           Vector3d(r*cos(theta+pi/2.0),
//				                                                    r*sin(theta+pi/2.0),
//				                                                    rr)));
//			}
//		}
//		for (auto it = ray_directions.begin(); it != ray_directions.end(); ++it) {
//			Vector3d ray_direction = (*it).first;
//			Vector3d ray_origin = (*it).second;
//			Ray ray(ray_origin, ray_direction);
//			std::list<Intersection> list_intersect = cone.hit(ray);
//			Intersection intersect = list_intersect.front();
//			std::pair<Vector3d,Vector3d> borders = intersect.get_path();
//			REQUIRE(list_intersect.size() == 0);
//		}
//	}
//
//	SECTION("One intersection at apex") {
//		std::list<std::pair<Vector3d,Vector3d>> ray_directions;
//		double theta;
//		for (int j = 0; j < 8; ++j) {
//			theta = j*pi/8.;
//			ray_directions.emplace_back(std::make_pair(Vector3d(cos(theta), sin(theta), 0.),
//			                                           Vector3d(cos(theta), sin(theta), 0.)));
//		}
//
//		for (auto it = ray_directions.begin(); it != ray_directions.end(); ++it) {
//
//			Vector3d ray_direction = (*it).first;
//			Vector3d ray_origin = (*it).second;
//			Ray ray(ray_origin, ray_direction);
//			std::list<Intersection> list_intersect = cone.hit(ray);
//			Intersection intersect = list_intersect.front();
//			std::pair<Vector3d, Vector3d> borders = intersect.get_path();
//			REQUIRE(list_intersect.size() == 1);
//			REQUIRE((borders.first - borders.second).isMuchSmallerThan(1E-06, 1E-06));
//			REQUIRE((borders.second - origin).isMuchSmallerThan(1E-06, 1E-06));
//		}
//	}
//
//	SECTION("One intersection at apex along border") {
//		std::list<Vector3d> ray_directions{Vector3d(1., 0., 1.),
//		                                   Vector3d(-1., 0., 1.),
//		                                   Vector3d(1., 0., -1.),
//		                                   Vector3d(-1., 0., 1.),
//		                                   Vector3d(0., 1., 1.),
//		                                   Vector3d(0., -1., 1.),
//		                                   Vector3d(0., 1., -1.),
//		                                   Vector3d(0., -1., -1.)};
//		Vector3d ray_origin(0., 0., 0.);
//		for (auto it = ray_directions.begin(); it != ray_directions.end(); ++it) {
//			Ray ray(ray_origin, *it);
//			std::list<Intersection> list_intersect = cone.hit(ray);
//			Intersection intersect = list_intersect.front();
//			std::pair<Vector3d,Vector3d> borders = intersect.get_path();
//			REQUIRE(list_intersect.size() == 1);
//			REQUIRE((borders.first + cone.big_scale()*ray.direction()).isMuchSmallerThan(1E-06, 1E-06));
//			REQUIRE((borders.second - cone.big_scale()*ray.direction()).isMuchSmallerThan(1E-06, 1E-06));
//		}
//
//	}
//}