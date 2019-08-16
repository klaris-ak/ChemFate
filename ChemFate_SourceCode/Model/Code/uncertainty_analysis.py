from collections import defaultdict
import copy
import numpy as np

chemProp_old = {'a': 1, 'b': 2}
uncertainty_params_old = {'a': 0.1, 'b': 0.2}

def sample_chem_prop(chemProp_old, uncertainty_params_old):
    chemProp_dict = defaultdict(list)
    chemProp = copy.deepcopy(chemProp_old)
    uncertaintyParams = copy.deepcopy(uncertainty_params_old)

    # calculate 1000 times triangular distribution, store in a params dict
    for _ in range(1000):
        for key in uncertaintyParams:
            if key in chemProp:
                value = chemProp[key]
                interval = uncertaintyParams[key]
                if value > 0 and interval > 0:
                    new_value = np.random.triangular(value * (1 - interval), value, value * (1 + interval))
                    chemProp_dict[key].append(new_value)

    # get the 1st, 2nd, ..., 99th percentile of the params
    percentile_dict = defaultdict(list)
    # for each parameter
    for k, v in chemProp_dict.iteritems():
        for i in range(100):
            percentile_dict[k].append(np.percentile(v, i))

    # use latin hypercube sampling, generate list of 100 dicts, each with same params
    sampling_list = [{} for _ in xrange(100)]
    print sampling_list

    total_size = 100
    for x in range(100):
        for k, v in percentile_dict.iteritems():
            index = np.random.randint(0, total_size)
            sampling_list[x][k] = v[index]
            v.pop(index)
        total_size -= 1

    return sampling_list


sample_chem_prop(chemProp_old, uncertainty_params_old)