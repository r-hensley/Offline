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
  //
  // Association lifecycle (implemented by KKStrawHit):
  //   updateState     -> save weight/times, advance generation, clear posterior;
  //   accepted fit    -> copy its posterior at the saved time for that generation;
  //   unbiased CA     -> use the pair only while activity, support and weight match.
  // An ordinary reference update clears the posterior without authorizing a new
  // one. The cache records this association; it performs no material transport.
  template<class KTRAJ> class KKStrawHitLocalReference {
    public:
      struct Cache {
        // A default entry can exist before any weight has been recorded.
        bool weightKnown = false;
        // These count weight updates for this object, not meta-iterations or
        // events. Zero posteriorGeneration denotes no accepted association.
        unsigned long long weightGeneration = 0;
        unsigned long long posteriorGeneration = 0;
        // Both CA clocks at weight creation, not the hit's latest biased CA.
        double particleTime = std::numeric_limits<double>::quiet_NaN();
        double sensorTime = std::numeric_limits<double>::quiet_NaN();
        // Exact measurement information (matrix AND vector) supplied to the fit.
        KinKal::Weights weight;
        // Owned post-fit piece at particleTime; never a pre-fit covariance cache.
        std::shared_ptr<KTRAJ> posterior;
        bool posteriorActive = false;
      };
      // Counts are completed calculations, including repeated calls for one hit.
      // Support/usability flags can overlap; they are not disjoint failure bins
      // or track-efficiency measures, and they do not test positive definiteness.
      struct Counts {
        std::size_t local = 0;
        std::size_t fallback = 0;
        std::size_t localOutsidePiece = 0;
        std::size_t localOutsideSensor = 0;
        std::size_t localUnusable = 0;
      };
    private:
      struct State {
        // Object addresses identify entries only within this scope. Constructors
        // call reset so allocator reuse cannot inherit a destroyed hit's entry.
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
      // Requires enabled(). A first lookup creates an ineligible default entry,
      // so inspection before the first weight update cannot authorize subtraction.
      static Cache& cache(void const* hit) { return current_->cache[hit]; }
      static void reset(void const* hit) {
        if(enabled()) current_->cache.erase(hit);
      }
      static bool sameWeight(KinKal::Weights const& a, KinKal::Weights const& b) {
        // This is an identity guard, not a numerical agreement test: the snapshot
        // must describe the exact stored information. A tolerance could silently
        // accept a different linearization or annealing scale. Check the vector
        // and every matrix entry, not merely variances or a projected DOCA error.
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
          // Reached only after the local inversion and CA return. An exception
          // before record() is neither a completed local call nor a fallback.
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
