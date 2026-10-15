#ifndef Mu2eKinKal_KKStrawHitLocalReference_hh
#define Mu2eKinKal_KKStrawHitLocalReference_hh

#include "KinKal/General/Weights.hh"
#include <cstddef>
#include <limits>
#include <memory>
#include <unordered_map>

namespace mu2e {
  // Experimental, scoped sidecar storage: it deliberately leaves the layout
  // and virtual interface of KKStrawHit unchanged. The owning module must keep
  // a Scope alive throughout fitting and any intended diagnostic sampling.
  template<class KTRAJ> class KKStrawHitLocalReference {
    public:
      struct Cache {
        bool weightKnown = false;
        unsigned long long weightGeneration = 0;
        unsigned long long posteriorGeneration = 0;
        double particleTime = std::numeric_limits<double>::quiet_NaN();
        double sensorTime = std::numeric_limits<double>::quiet_NaN();
        KinKal::Weights weight;
        std::shared_ptr<KTRAJ> posterior;
        bool posteriorActive = false;
      };
      struct Counts {
        std::size_t local = 0;
        std::size_t fallback = 0;
        std::size_t localOutsidePiece = 0;
        std::size_t localOutsideSensor = 0;
        std::size_t localUnusable = 0;
      };
    private:
      struct State {
        std::unordered_map<void const*,Cache> cache;
        Counts counts;
      };
      inline static thread_local State* current_ = nullptr;
    public:
      // Each call starts empty, including a nested call. Restoring the previous
      // scope on destruction prevents exceptions or address reuse across events
      // from attaching a posterior to a different hit. Disabled scopes allocate
      // no entries and explicitly suspend any enclosing diagnostic scope.
      class Scope {
        State state_;
        State* previous_;
        public:
          explicit Scope(bool enabled) : previous_(current_) {
            current_ = enabled ? &state_ : nullptr;
          }
          ~Scope() { current_ = previous_; }
          Scope(Scope const&) = delete;
          Scope& operator=(Scope const&) = delete;
          Counts const& counts() const { return state_.counts; }
      };
      static bool enabled() { return current_ != nullptr; }
      static Cache& cache(void const* hit) { return current_->cache[hit]; }
      static void reset(void const* hit) {
        if(enabled()) current_->cache.erase(hit);
      }
      static bool sameWeight(KinKal::Weights const& a, KinKal::Weights const& b) {
        for(unsigned i=0; i<KinKal::NParams(); ++i) {
          if(a.weightVec()[i] != b.weightVec()[i]) return false;
          for(unsigned j=0; j<KinKal::NParams(); ++j)
            if(a.weightMat()(i,j) != b.weightMat()(i,j)) return false;
        }
        return true;
      }
      static void record(bool local, bool inPiece, bool inSensor, bool usable) {
        if(!enabled()) return;
        auto& counts = current_->counts;
        if(local) {
          ++counts.local;
          if(!inPiece) ++counts.localOutsidePiece;
          if(!inSensor) ++counts.localOutsideSensor;
          if(!usable) ++counts.localUnusable;
        } else {
          ++counts.fallback;
        }
      }
  };
}
#endif
