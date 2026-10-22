// Standalone regression for the six-float BkgANNSHU inference contract.
// Uses only synthetic coefficients; no trained model or event data is needed.
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>
#include <sys/mman.h>
#include <unistd.h>

#include "Offline/Mu2eKinKal/inc/TrainBkg.hxx"

namespace {
  class TemporaryWeights {
    public:
      TemporaryWeights() {
        char pattern[] = "/tmp/trainbkg-single-hit-XXXXXX";
        int fd = mkstemp(pattern);
        if (fd < 0) throw std::runtime_error("mkstemp failed");
        path_ = pattern;
        close(fd);
      }
      ~TemporaryWeights() { std::remove(path_.c_str()); }
      std::string const& path() const { return path_; }
    private:
      std::string path_;
  };

  class GuardedInput {
    public:
      GuardedInput() {
        long page = sysconf(_SC_PAGESIZE);
        if (page <= 0) throw std::runtime_error("sysconf failed");
        size_ = 2 * static_cast<std::size_t>(page);
        region_ = mmap(nullptr, size_, PROT_READ | PROT_WRITE,
                       MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
        if (region_ == MAP_FAILED) throw std::runtime_error("mmap failed");
        auto boundary = static_cast<char*>(region_) + page;
        if (mprotect(boundary, page, PROT_NONE) != 0) {
          munmap(region_, size_);
          throw std::runtime_error("mprotect failed");
        }
        // The seventh float is on an inaccessible page. This also detects
        // reads inside BLAS libraries built without address sanitizers.
        input_ = reinterpret_cast<float*>(boundary) - 6;
      }
      ~GuardedInput() { munmap(region_, size_); }
      float* data() const { return input_; }
    private:
      void* region_;
      std::size_t size_;
      float* input_;
  };

  void writeTensor(std::ostream& out, std::string const& name,
                   std::vector<float> const& values) {
    out << name << ' ' << values.size() << '\n';
    for (float value : values) out << value << ' ';
    out << '\n';
  }

  void writeWeights(std::string const& path) {
    std::ofstream out(path);
    if (!out) throw std::runtime_error("opening synthetic weights failed");
    // Split every feature into positive/negative ReLU channels, pass them
    // through two diagonal layers with negative biases, then recombine.
    // The resulting channels are max(+/- x - 0.625, 0).
    std::vector<float> output(12), identity(144, 0.0f), split(72, 0.0f);
    for (int i = 0; i < 6; ++i) {
      output[2*i] = (i + 1) / 16.0f;
      output[2*i + 1] = -output[2*i];
      split[i*12 + 2*i] = 1.0f;
      split[i*12 + 2*i + 1] = -1.0f;
    }
    for (int i = 0; i < 12; ++i) identity[i*12 + i] = 1.0f;
    writeTensor(out, "tensor_dense3bias0", {0.125f});
    writeTensor(out, "tensor_dense3kernel0", output);
    writeTensor(out, "tensor_dense2bias0", std::vector<float>(12, -0.5f));
    writeTensor(out, "tensor_dense2kernel0", identity);
    writeTensor(out, "tensor_dense1bias0", std::vector<float>(12, -0.25f));
    writeTensor(out, "tensor_dense1kernel0", identity);
    writeTensor(out, "tensor_densebias0", std::vector<float>(12, 0.125f));
    writeTensor(out, "tensor_densekernel0", split);
    out.close();
    if (!out) throw std::runtime_error("writing synthetic weights failed");
  }
}

int main() {
  try {
    TemporaryWeights weights;
    writeWeights(weights.path());
    TMVA_SOFIE_TrainBkg::Session session(weights.path());
    // The session has loaded the coefficients. Remove the file before the
    // guarded call so even a regression that faults leaves no temporary asset.
    if (std::remove(weights.path().c_str()) != 0)
      throw std::runtime_error("removing synthetic weights failed");
    GuardedInput input;
    std::vector<std::array<float, 6>> cases{
      {{0, 0, 0, 0, 0, 0}},
      {{-1, 2, -3, 4, -5, 6}},
      {{1, 1, 1, 1, 1, 1}},
      {{-0.5f, 0.5f, -0.5f, 0.5f, -0.5f, 0.5f}},
      {{-0.625f, 0.625f, -0.125f, 0.125f, -0.25f, 0.25f}},
      {{100, 100, 100, 100, 100, 100}},
      {{-100, -100, -100, -100, -100, -100}}
    };
    for (int feature = 0; feature < 6; ++feature) {
      for (float value : {-2.0f, 2.0f}) {
        std::array<float, 6> impulse{};
        impulse[feature] = value;
        cases.push_back(impulse);
      }
    }
    cases.push_back(cases.front()); // Reset after nonzero inputs.
    float maxError = 0;
    // Reuse the same session to exercise all of its intermediate buffers.
    for (auto const& values : cases) {
      std::copy(values.begin(), values.end(), input.data());
      auto result = session.infer(input.data());
      if (result.size() != 1) throw std::runtime_error("expected one score");
      float logit = 0.125f;
      for (int i = 0; i < 6; ++i)
        logit += std::copysign(std::max(std::abs(values[i]) - 0.625f, 0.0f), values[i])
          * (i + 1) / 16.0f;
      float expected = 1.0f / (1.0f + std::exp(-logit));
      if (!std::isfinite(result[0]) || std::abs(result[0] - expected) > 1e-6f)
        throw std::runtime_error("synthetic score mismatch");
      if (!std::equal(values.begin(), values.end(), input.data()))
        throw std::runtime_error("input was modified");
      maxError = std::max(maxError, std::abs(result[0] - expected));
    }
    std::cout << "TrainBkg single-hit guarded inference passed: " << cases.size()
      << " cases, max absolute error " << maxError << '\n';
  } catch (std::exception const& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
  return 0;
}
