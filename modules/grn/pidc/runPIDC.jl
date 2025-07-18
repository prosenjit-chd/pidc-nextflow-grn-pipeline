using NetworkInference
using LightGraphs

algorithm = PIDCNetworkInference()

dataset_name = string(ARGS[1])  # Input TSV file

@time genes = get_nodes(dataset_name);
@time network = InferredNetwork(algorithm, genes);

write_network_file(string(ARGS[2]), network);
