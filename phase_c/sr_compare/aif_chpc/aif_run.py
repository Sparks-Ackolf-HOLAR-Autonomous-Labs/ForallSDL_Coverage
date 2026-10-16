import json, os, time
import aifeynman

os.chdir(os.path.dirname(os.path.abspath(__file__)))
out = {}
for task in ["literature", "md", "combined"]:
    t0 = time.time()
    try:
        aifeynman.run_aifeynman("data/", f"{task}.txt", 30,
                                "14ops.txt", polyfit_deg=3, NN_epochs=400)
        sol = f"results/solution_{task}.txt"
        lines = open(sol).read().strip().splitlines() if os.path.exists(sol) else []
        out[task] = {"seconds": round(time.time() - t0, 1),
                     "solution_lines": lines[-8:]}
    except Exception as e:
        out[task] = {"seconds": round(time.time() - t0, 1), "error": repr(e)}
    with open("aif_out.json", "w") as f:
        json.dump(out, f, indent=1)
print(json.dumps(out, indent=1))
