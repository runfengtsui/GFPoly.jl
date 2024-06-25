# GFPoly.jl

判断 Galois 域上的多项式是否为不可约多项式或者本原多项式. 同时暴力查找列出所有 $n$ 次不可约多项式或本原多项式.

## 不可约多项式判定定理

给定有限域为 $\mathbb{F}_p$, $\mathbb{F}_p[x]$ 为 $\mathbb{F}_p$ 上的一元多项式环, $f(x)\in\mathbb{F}_p[x]$ 为 $n$ 次多项式, 若对任意的 $l$, $l\in\mathbb{Z}^+$, $1\le l\le\left[\dfrac{n}{2}\right]$, 都有 $gcd \left(x^{p^l}-x,f(x)\right)=1$, 那么 $f(x)$ 为 $n$ 次不可约多项式.

## 本原多项式判别定理

给定有限域 $\mathbb{F}_p$, $\mathbb{F}_p[x]$ 为 $\mathbb{F}_p$ 上的一元多项式环, $f(x)\in\mathbb{F}_p[x]$ 为 $n$ 次($n\ge 2$)多项式, 根据整数的唯一分解定理得: $N=p^n-1=p_1^{l_1}p_2^{l_2}\cdots p_w^{l_w}$, 若满足

(1) 任一正整数 $l \left(1\le l\le\left[\dfrac{n}{2}\right]\right)$, 都有 $gcd \left(x^{p^l}-x,f(x)\right)=1$;

(2) 任一正整数 $i(1\le i\le w)$, 都有 $f(x) \nmid x^{\frac{p^n-1}{p_i}}-1$;

那么 $f(x)$ 为 $n$ 次本原多项式.

## 参考资料

* [GaloisFields.jl](https://github.com/tkluck/GaloisFields.jl)
* [AbstractAlgebra.jl](https://github.com/Nemocas/AbstractAlgebra.jl)
* [GaussianPrimeField](https://github.com/Bocconi-StatCSPhD-CS1-2023/GaussianPrimeField)
* 裴定一, 祝跃飞. 算法数论 | 2版[M]. 北京：科学出版社, 2015
* 杨福祥. 有限域上本原多项式的研究[D]. 上海交通大学, 2011.

