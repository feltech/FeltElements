//
// Created by dave on 04/11/2019.
//

#ifndef FELTELEMENTS_MESHFILE_HPP
#define FELTELEMENTS_MESHFILE_HPP

#include <memory>
#include <string>

class tetgenio;

namespace FeltElements
{
class MeshFile
{
  public:
	explicit MeshFile(std::string file_name);
	~MeshFile();

	[[nodiscard]] std::size_t num_simplexes() const;
	[[nodiscard]] std::size_t num_corners() const;
	[[nodiscard]] std::size_t num_trifaces() const;
	[[nodiscard]] double * points() const;
	[[nodiscard]] double * displacements() const;
	[[nodiscard]] int * corners() const;

	[[nodiscard]] double * tet_corner_point(std::size_t tet_idx, std::size_t corner_idx) const;
	[[nodiscard]] double *
	tet_corner_displacement(std::size_t tet_idx, std::size_t corner_idx) const;

  private:
	std::string m_file_name;
	std::unique_ptr<tetgenio> m_io;
};
} // namespace FeltElements

#endif // FELTELEMENTS_MESHFILE_HPP
