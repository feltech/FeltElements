//
// Created by dave on 04/11/2019.
//

#ifndef CPPLUA_MESHFILE_HPP
#define CPPLUA_MESHFILE_HPP

#include <memory>
#include <string>

class tetgenio;

namespace FeltElements
{
class MeshFile
{
  public:
	explicit MeshFile(std::string file_name);
	virtual ~MeshFile();

	[[nodiscard]] std::size_t num_simplexes() const;

  private:
	std::string m_file_name;
	std::unique_ptr<tetgenio> m_io;
};
}

#endif // CPPLUA_MESHFILE_HPP
