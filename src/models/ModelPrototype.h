#pragma once

#include <Model.h>

#include <memory>

namespace miho {

enum class ModelType { geo, abc };

class ModelPrototype : public Model {
 public:
  virtual ~ModelPrototype() noexcept = default;
  /// The prototype design pattern
  virtual std::unique_ptr<ModelPrototype> clone() const = 0;
  /// external interface
  virtual void set_sigma(const std::vector<double>& sigma) = 0;
};

}  // namespace miho
