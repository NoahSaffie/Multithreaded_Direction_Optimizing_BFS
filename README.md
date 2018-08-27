# Multithread_BFS

Notes: 
Some optimization is loss due to the need for omp critical pragma when updating neighbors in Top Down BFS. Without the critical pragma, speed is greatly increased however output is incorrect due to the possibility of a neighbor being visited more than once cuncurrently and marked as so.
