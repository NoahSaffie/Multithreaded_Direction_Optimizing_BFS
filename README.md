# Multithread_Asynchronous_BFS

Notes: 
Some optimization is loss due to the need for omp critical pragma when updating neighbors in Top Down BFS. Without the critical pragma, speed is greatly increased however output is incorrect due to the possibility of a neighbor being visited more than once cuncurrently and marked as so. Most importantly, because of the condition of the main loop, an incorrect reporting/recording of node amounts can cause the program to abort early (since it thinks it found all the correct nodes).
