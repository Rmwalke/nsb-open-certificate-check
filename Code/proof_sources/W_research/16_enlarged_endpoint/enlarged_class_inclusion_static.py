"""Fresh exact class-inclusion and witness checks.

The endpoint certificate is read as data; no previous program is imported
or executed. The comparisons concern the two named earlier analytic classes.
"""
from fractions import Fraction as F
from pathlib import Path
import hashlib
import json


def main():
    path = Path(__file__).with_name("enlarged_endpoint_gate_result.json")
    raw = path.read_bytes()
    data = json.loads(raw)
    checks = []

    def check(label, condition):
        assert condition, label
        checks.append(label)

    check("endpoint certificate passed",
          data["status"] == "PASS_ENLARGED_ENDPOINT_REGULARITY_CLASSES")
    passed = {(F(c["amplitude"]), F(c["chi"]))
              for c in data["finite_candidate_classes"] if c["passed"]}
    old = [("old_A", F(11,10), F(1,100)),
           ("old_B", F(1), F(1,20))]
    new_amplitude = (F(8,7), F(1,100))
    new_radius = (F(1), F(1,10))
    check("larger amplitude class certified", new_amplitude in passed)
    check("larger radius class certified", new_radius in passed)
    check("old A contained in new amplitude class",
          old[0][1] <= new_amplitude[0] and old[0][2] <= new_amplitude[1])
    check("old B contained in new radius class",
          old[1][1] <= new_radius[0] and old[1][2] <= new_radius[1])

    # W has six unit-frequency Fourier vectors, each of Euclidean norm 1/2.
    W_chi = 6*F(1,2)
    W_E = 6*F(1,2)**2
    check("W critical norm exactly three", W_chi == 3)
    check("W energy and enstrophy exactly three halves", W_E == F(3,2))

    # h=(0,(1/5)sin(2x),0): two modes at |k|=2, each of norm 1/10.
    h_amplitude = F(1,5)
    h_coefficient = h_amplitude/2
    h_chi = 2*h_coefficient/2
    h_E = 2*h_coefficient*h_coefficient
    h_X = 4*h_E
    check("higher-shell h critical norm one tenth", h_chi == F(1,10))
    check("higher-shell h energy one fiftieth", h_E == F(1,50))
    check("higher-shell h enstrophy two twenty-fifths", h_X == F(2,25))
    check("higher-shell witness admitted", (F(1), h_chi) in passed)

    witnesses = []
    for name, amplitude, extra_chi, energy, enstrophy, field in [
        ("higher_shell", F(1), h_chi, W_E+h_E, W_E+h_X,
         "(sin z, sin x+(1/5)sin(2x), sin y)"),
        ("larger_amplitude", F(8,7), F(0),
         F(8,7)**2*W_E, F(8,7)**2*W_E,
         "(8/7)*(sin z,sin x,sin y)"),
    ]:
        # Disjoint supports give chi(u0-bW)=3*|amplitude-b|+extra_chi.
        # Its minimum over |b|<=A is 3*max(amplitude-A,0)+extra_chi.
        exclusions = []
        for label, maximum, radius in old:
            distance = W_chi*max(amplitude-maximum,F(0))+extra_chi
            check(f"{name} outside every decomposition in {label}",
                  distance > radius)
            exclusions.append({
                "old_class": label, "old_amplitude_cap": str(maximum),
                "old_critical_radius": str(radius),
                "minimum_critical_distance_to_old_segment": str(distance),
                "strict_gap": str(distance-radius),
            })
        witnesses.append({
            "name": name, "initial_field": field,
            "smooth_real_mean_zero_divergence_free": True,
            "divergence_reason": "Each component is independent of its matching coordinate.",
            "E0": str(energy), "X0": str(enstrophy),
            "distance_identity": f"A_-1(u0-bW)=3*abs({amplitude}-b)+{extra_chi}",
            "outside_both_old_classes_under_any_allowed_decomposition": True,
            "exclusions": exclusions,
        })

    result = {
        "status": "PASS_STRICT_ENLARGEMENT_OF_TWO_PRIOR_ANALYTIC_CLASSES",
        "source_json": path.name,
        "source_sha256": hashlib.sha256(raw).hexdigest(),
        "old_classes": [{"name": n, "amplitude": str(a), "chi": str(r)}
                        for n,a,r in old],
        "new_amplitude_class": {"amplitude": "8/7", "chi": "1/100"},
        "new_critical_radius_class": {"amplitude": "1", "chi": "1/10"},
        "critical_radius_enlargement_factor": "2",
        "amplitude_endpoint_increment": str(F(8,7)-F(11,10)),
        "old_union_is_strict_subset_of_new_union": True,
        "witnesses": witnesses,
        "assertion_count": len(checks),
        "checks": checks,
        "scope": "Strict comparison with the two explicitly named earlier analytic classes, not a novelty claim or a claim about every conceivable prior sufficient criterion.",
    }
    out = path.with_name("enlarged_class_inclusion_static_result.json")
    out.write_text(json.dumps(result,indent=2)+"\n",encoding="utf-8")
    print(json.dumps({"status":result["status"],"assertions":len(checks),
                      "witnesses":[w["initial_field"] for w in witnesses]},
                     indent=2))


if __name__ == "__main__":
    main()
