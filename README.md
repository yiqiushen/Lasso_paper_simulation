# Lasso_paper_simulation
R code using vanilla CVXR to reproduce simulation results in contaminated lasso paper (see PDF).

# Why not `glmnet`?
In our paper, we are exploring the performance of Lasso under corrupted response and correlated design with a fixed penalization $\lambda$. On the other hand, `` `glmnet` fits the entire solution path for the lasso or elastic net problems efficiently with various techniques such as using warm starts and strong rules. Those advantages will disappear if the $\lambda$ sequence is forced to be only one value'' (Ref: https://glmnet.stanford.edu/articles/glmnet.html#appendix-2-comparison-with-other-packages). In fact, when we run experiments using `glmnet`, even when $\lambda$ is set to 0, a situation where Lasso with high probability guarantees exact recovery, the prediction risk is bounded away from 0. 


