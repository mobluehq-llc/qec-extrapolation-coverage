# F2 repro: rotated surface code memory-Z, stim+pymatching via sinter. Compare to 0.1*(100p)^((d+1)/2) (Fowler/Gidney-Ekera 2019 scaling).
import sys, json, stim, sinter
def main():
    model = sys.argv[1]          # 'full' (uniform circuit-level depol) or 'dataonly' (negative control)
    ps = [1e-3, 2e-3, 3e-3, 5e-3]
    tasks = []
    for d in (3, 5, 7):
        for p in ps:
            if model == 'full':
                c = stim.Circuit.generated("surface_code:rotated_memory_z", distance=d, rounds=d,
                    after_clifford_depolarization=p, before_round_data_depolarization=p,
                    before_measure_flip_probability=p, after_reset_flip_probability=p)
            else:  # data-qubit depolarization only, perfect gates/measurement
                c = stim.Circuit.generated("surface_code:rotated_memory_z", distance=d, rounds=d,
                    before_round_data_depolarization=p)
            tasks.append(sinter.Task(circuit=c, decoder='pymatching', json_metadata={'d': d, 'p': p, 'r': d, 'model': model}))
    stats = sinter.collect(num_workers=4, tasks=tasks, max_shots=int(sys.argv[2]), max_errors=int(sys.argv[3]), print_progress=False)
    out = [dict(d=s.json_metadata['d'], p=s.json_metadata['p'], r=s.json_metadata['r'], shots=s.shots, errors=s.errors, model=model) for s in stats]
    json.dump(out, open(f'repro_{model}.json', 'w'), indent=1)
    print('done', len(out))

if __name__=='__main__':
    main()
