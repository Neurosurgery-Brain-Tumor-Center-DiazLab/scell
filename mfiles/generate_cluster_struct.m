function cluster = generate_cluster_struct(score)

cluster = struct;
cluster.my_list_1 = randi(2, size(score,1), 1);
cluster.my_list_2 = randi(4, size(score,1), 1);
cluster.my_list_3 = randi(6, size(score,1), 1);