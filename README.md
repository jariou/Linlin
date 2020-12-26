[[_TOC_]]

# Optimal Divid of The Compound Poisson Process under a Stochastic Interest Rate

 Linlin Tiana† Lihua Baib‡ Junyi Guoc

 a. School of Mathematical Sciences, Nankai University, Tianjin 300071, China.
 
b. School of Mathematics Sciences, Nankai University, Tianjin 300071, China.

c. School of Mathematics Sciences, Nankai University, Tianjin 300071, China.

**Abstract** In this paper we assume the insurance wealth process is driven by the compound poisson process. The discounting factor is modelled as a geometric[1] Brownian motion at first and then as an exponential function of an integrated Ornstein-Uhlenbeck process. The objective is to maximize the cumulated value of expected discounted dividends up to the ruin time[2]. We give an explicit expression of the value function and the optimal strategy in the first case. For the second case, we explore[3] some properties of value function[4]. Since we can not find the explicit expression of the value function[5] in the second case, we use the notion of viscosity solution and prove that the value function is the viscosity solution of the corresponding HJB equation.

**Keywords:** Hamilton-Jacobi-Bellman equation, 
Vasicek model, 
Geometric Brownian motion, 
Interest rate, 
Viscosity solution, 
Optimal dividends

## Introduction
The optimal dividend problem has been discussed for a long time in the literature. In 1957 De Finetti [8] proposed that the insurance company should allow cash leakages and measure their performance during its life time instead of only focussing on ruin probability.

These cash leakages can be interpreted as dividends. In the setting of constant interest rate, Asmussen and Taksar (1997) [4] solved the optimal dividend problem for the special case of Brownian motion. They found out that the optimal strategy is a constant barrier strategy in the case of unbounded dividend and a so-called threshold strategy in the case of restricted dividend rates. In the case of surplus process following the compound poisson process, Gerber and Shiu (2006) [11] showed that the optimal strategy is a threshold strategy when claim size are exponentially distributed for restricted dividend rates. For more general case of claim size distribution, Azcue and Muler (2005) [3] studied the optimal reinsurance and dividend policy in the frame work[6] of Cram´er-Lundberge[7] model using viscosity solution. Later[8] in 2012, they found the optimal dividend payment policy which maximizes the cumulative expected discounted dividend pay-outs up to ruin time with a ceiling on the dividend rates. In the setting of constant interest rate the optimal dividend problem has been studied quite well under various general reserve models, see e.g., [1, 12, 15]. We omit listing the existing literature and refer to a survey on the dividend problems by Albrecher and Thonhauser (2009) [2] and references therein.

Interest rate forms a key component of the financial market, influencing on the firm’s cost and profit[9]. As it changes over time, it’s more reasonable to assume that the interest rates is a function[10] of time instead of a deterministic constant. The changes of interest rate reflects the fluctuations of the monetary market. Eisenberg [9] published a paper about optimal dividends in the setting of surplus follows[11] a drifted Brownian motion in 2015. The interest rate in this model also follows a Brownian motion with drift. They found an explicit expression for the value function of the optimal strategy for both restricted and unrestricted dividends. In section 2 of our paper, we model the surplus process as a compound poisson process and the claim size follows exponential distribution[12] and the interest rate follows a Brownian motion with drift. We explored dividend maximization problem[13] under the Dothan model and find[14], similar to the case of deterministic interest rate, that the optimal strategy does not change (compared to the Gerber-Shiu case) in its form, but in the parameters[15]. In Section 3, we consider Vasicek model, which the short rate[16] is defined as an Ornstein-Uhlenbeck process. Here, the situation changed completely. It is not that easy to calculate the return function of the corresponding strategy. We explore the continuity of value function but unfortunately we can not prove more regularities[17] about the value function. It’s natural to consider the problem in the framework of viscosity solutions. 

We recall that the notion of viscosity solutions was introduced by Crandall and Lions [7] for the first order equations and Lions [13, 14] for the second order equations. The notion of viscosity solution of integro-differential equations was pursued by Soner(1986) [16]. The viscosity solution concept of fully nonlinear partial differential equations has been proved to be extremely useful for the[18] control theory due to the fact that it do[19] not need the differentiability of the value function. It merely requires the[20] continuity of value function to define the viscosity solution. We refer to the user’s guide of Crandall, Ishii and Lions [6] for an overview of the theory of viscosity solutions and their applications. Using the notion of viscosity solution we characterize the value function as a viscosity solution of the corresponding Hamilton-Jacobi-Bellman equation.

## Geometric Brownian motion as a discounting factor

### Problem Formulation
Consider an insurance company whose surplus is given by a compound poisson process with a constant premium rate. The insurance company is allowed to pay out dividends, where the accumulated dividends until $t$ are given by $L_t$. The surplus at time $t$ is described as:
$$
X^L = x + ct 
$$
and
$$
- \sum^{N(t)}_{i=k} Y_k - L_t,
$$
where $x \ge 0$ is the initial surplus, the constant $c \gt 0$ is the premium rate, $N(t)$ is the poisson process representing the frequency of the incoming claims, whereas 
$\left{ Y_i \right}_{ i = 1 }^\infty$ 
is a sequence of random variables representing size of the incoming claims. We assume that Yi are independent, identically distributed (i.i.d.) with a common exponential intensity function $p(y) = βe^{-βy}$, where $β$ is a positive constant. Let $r_t = r + mt + δB_t$, where $m$ and $δ$ are two constants and $m > \frac{δ^2}{2}$. $B_t$ is a standard Brownian motion independent of $X^L_t$. As a risk measure, we consider that dividends are discounted by the geometric Brownian motion 

$$
exp\{-r - mt - δBt\}.
$$ 

Let 
$$\left(\Omega, F,(F_t)_{t\ge 0}, P\right)$$ 
be a complete probability space generated by $\{X_t, r_t\}$. Here we only allow the restricted dividend, which means, the cumulative dividend up to time $t$ is given by 
$L_t = \int_0^t l_s ds$, 
with $l_s \in [0, M]$ for some $0 < M < c$. We say that a strategy $L$ is admissible if it is predictable, nondecreasing, cadlag and it verifies $X^L_t \ge 0$ up to the ruin time. Denotes $U_{ad}[r, x]$ the set of all admissible strategies. 

Denoting by $τ^L$ L the ruin time of the surplus process under some admissible strategy L = {ls}, we define the return function corresponding to $L$ to be 
$$
J ^L(r, x) = E\left[\int^{ τ^L}_0 e^{-r-ms-δB_s}l_sds\right]. \ \ \ \ \ \ \ (2.1)
$$

And our objective is to find out the optimal strategy $L$ such that
$$
V (r, x) = {\sup J^L(r, x) \atop L\in Uad[r,x]}.\ \ \ \ \ \ \ \ \ \ \ \ \ \ \  (2.2)
$$
Here we note that for any strategy $L$,
$$
J^L(r, x) \le 
				E\left[\int^{τ^L}_0 e^{-r-ms-δB_s} c_s ds\right] 
				\le 
				E\left[\int^{τ^L}_0 e^{-r-ms-δB_s} M ds\right] 
				= 
				\frac{Me^{-r}} {m -\frac{δ^2}{2}}.
$$

This means $V (r, x)$ is bounded. 

According to Fleming and Soner (2006) [9], the HJB equation corresponding to the problem is
$$
\left[mVr + \frac{δ^2}{2}V_{rr} + 
cV_x - λV \right](r, x) + 
λ\int ^x_0 V (r, x - y)βe^{-βy}dy + 
{max \atop l \in[0,M]}(e^{-r} - V_x(r, x))l = 
0. \ \ \ \ \ \ \ \ \ \ \ \ (2.3)
$$

### Verification of Optimality

**Proposition 2.1** Assume there exists a function 
$v(r, x) \in C^{2,1}(R \times [0, +\infty))$ 
satisfies the HJB equation (2.3). For any dividend strategy {L_t}, we claim that 

$$
v(r, x) \ge E\left[\int^{τ^L}_0 e^{-r_t} dl(t)|r_0 = r, x_0 = x\right], (2.4)
$$

where $τ^L$ denotes the ruin time of insurance company with the dividend strategy $L_t$. Further, if there exists a strategy $L^{\star}_t$ such that $v(r, x) = E\left[\int^{τ^L}_0 e^{-r_t}dL^{\star}_t|r_0 = r, x_0 = x\right]$, then strategy $L^{\star}_t$ must be optimal.

### Solve V(r,x;b)
For a constant threshold b ≥ 0, Let Xbt denotes the surplus process with threshold strategy L b and τ b denotes the corresponding ruin time of the insurance company. We can see V (r, x; b) = E[Z τb0e−r−ms−δBs dls] = e−rE[Z τb0e−ms−δBs dls]. (2.7)

Since τ b is independent of rt, we can denote F(x; b) := E[R τb0e−ms−δBs dls], which implies Vx(r, x; b) = e−rF0(x; b), Vxx(r, x; b) = e
−rF 00(x; b), Vr(r, x; b) = −e −rF(x; b). Substituting e−rF(x; b) into (2.5), we can obtain that F(x; b) satisfies the following equation
(−m +δ 2 2− λ)F(x; b) + cF0(x; b) + λZ x0 F(x − y; b)βe−βydy = 0, 0 < x < b, (2.8)
(−m +δ22− λ)F(x; b) + cF0 (x; b) + λ Z x 0F(x − y; b)βe−βydy + M − MF0
(x; b) = 0, x > b.(2.9)

Combining (2.3), (2.5) and V (r, x; b) = e −rF(x; b), we can see that if the optimal strategy is the threshold strategy with b∗ > 0, we must have F 0 (x; b ∗ ) > 1, for x < b∗ , (2.10)
F 0 (x; b ∗ ) < 1, for x > b∗. (2.11)

Following the similar way of Gerber and Shiu [11] in 2006, we can obtain the explicit expression of F(x; b). Equation (2.8) can be rewritten as cF00(x; b) + (βc − λ − (m − 1 2 δ 2 ))F 0 (x; b) − β(m −1
2δ2)F(x; b) = 0, (2.12)
with a general solution of the form F(x; b) = C0e rx + C1e sx, (2.13)
where r > 0, s < 0 are the roots of the characteristic equation:
cξ2 + (βc − λ − (m − 1 2 δ 2 ))ξ − β(m − 1 2 δ 2 ) = 0, (2.14)
from the fact that c(−β) 2 + (βc − λ − (m − 1 2 δ 2 ))(−β) − β(m −
1 2 δ 2 ) > 0, we can see that s + β is positive. Substituting (2.13) into equation (2.8) and equating the coefficient e−βx with 0, we have
λβ( C0 r + β + C1 s + β ) = 0. We can rewrite F(x; b) = γ[(r + β)e
rx − (s + β)e sx], 0 ≤ x ≤ b, (2.15)
where γ does not depend on x. Similarly, we can rewrite (2.9) as
(c − M)F 00(x; b) + (β(c − M) − λ − (m − 1 2 δ 2 ))F 0 (x; b) − β(m −
1 2 δ 2 )F(x; b) + βM = 0, (2.16)
5
and deduce that (2.9) has a solution of the form
F(x; b) = M m − 1 2 δ 2 + Deux, x > b, (2.17)
where u is the negative solution of
(c − M)ξ 2 + [β(c − M) − λ − (m − 1 2 δ 2 )]ξ − β(m − 1 2  δ 2 ) = 0.
Here we can notice F(b; b) = M m − δ 2 2 + Deub . (2.18)

By the continuity of function, F(b−; b) = F(b+; b), we have γ[(r + β)e
rb − (s + β)e sb] = M m − 1  δ 2 + Deub . (2.19)
We substitute (2.15) and (2.17) into (2.9), set the coefficient e
−βx to 0 and then cancel the
factor βeβb to obtain γ(e rb − e sb) − M β(m − 1 2 δ 2 ) − Deub β + u
= 0. (2.20)
Combining (2.19) and (2.20), we obtain that γ = −   β M m − 1 2 δ 2 1
(r − u)e rb − (s − u)e sb . (2.21) Substituting (2.21) into (2.15), we obtain F(x; b) = −u β M m − 1 2  2 (β + r)e rx − (β + s)e
sx (r − u)e rb − (s − u)e sb , 0 ≤ x ≤ b. (2.22) Here we notice that
F(0; 0) = −u β M  m − 1 2 δ 2 . (2.23)  Combining (2.18) with (2.17), we obtain F(x; b) = M m − 1  2 δ 2 [1 − e u(x−b) ] + F(b; b)e u(x−b)
, x ≥ b. (2.24) Then we get F 0 (x; b) = −u[ M m − 1 2 δ 2 − F(b; b)]eu(x−b) , x ≥ b. (2.25) If b ∗ = 0, it means for all x ≥ 0, F
0 (x; 0) ≤ 1 or equivalently, if (−u)[ M m − 1 2 δ 2 − F(0; 0)] = (−u) M m − 1 2 δ 2 (1 + u β ) ≤ 1, then b ∗ = 0. 6 If (−u) M m− 1 2 δ 2 (1 + u β ) > 1, we can use the fact that F0 (b ∗−; b ∗ ) = 1 or F 0 (b ∗+; b ∗ ) = 1 to obtain the closed-form expression of the optimal threshold b ∗ . It turns out that b ∗ = 1 r − s ln(s 2 − us r 2 − ur ).
Using the fact F 0 (b ∗ ; b ∗ ) = 1 and (2.25) we can obtain F(b ∗ ; b
∗) = M m − 1 2 δ 2 + 1 u. (2.26)
Combining (2.26) and (2.24), we obtain F(x; b ∗ ) = M m − 1 2 δ 2 + 1 u e u(x−b ∗) , x > b∗ . (2.27)
As a summary, we give out the following proposition.
Proposition 2.2 The value function is given by V (r, x) = e
−rF(x; b ∗ ). If (−u) M m− 1 δ2(1 + u   β) ≤ 1,F(x; b∗) = M m −1
2 δ 2 (1 − e ux(1 + u β )), and b ∗ = 0. If (−u) M m− 1 2 δ 2 (1 + u
β ) > 1, F(x; b ∗ ) =    −u β M m− 1 2 δ 2 (β+r)e rx−(β+s)e sx (r−u)e rb∗−(s−u)e sb∗ , x ≤ b ∗  M m− 1 2 δ 2 + 1 u e u(x−b ∗) , x > b∗
.
(2.28)
and b ∗ = 1 r−s log( s 2−us r 2−ur). The optimal dividend strategy L
∗ = {l ∗ s} is l ∗ s = M1{XL∗ t >b∗} , here 1{XL∗ t >b∗} is the indicator function.

proof We only need to show that F 0 (x; b ∗ ) > 1 for x < b∗ and F
0 (x; b ∗ ) < 1 for x > b∗ . From (2.28), it’s very easy to verify that F0(x; b ∗) = eu(x−b∗) < 1 for x > b∗
.
To show (2.10) holds, we only need to show that F
00(x; b ∗ ) < 0 for 0 ≤ x < b∗. Differentiating (2.15), we obtain
F 00(x; b ∗ ) = γ[r 2 (r + β)e rx − s 2 (s + β)e sx], 0 < x < b∗.
It’s easy to deduce that F00(x; b∗) attains it’s maximum at x = b
∗. Consequently, we onlyneed to show F00(b∗−; b∗) ≤ 0. From (2.12), we see cF00(b ∗−; b ∗ ) = −(βc − λ − (m − 1 2 δ 2 )) + β(m − 1 2 δ 2 )F(b
∗; b∗),7 where Lt denotes the cumulated dividend process, c > 0 denotes the premium rate as in section before. Here we only allow the restricted dividend, which can be written as Lt = R t 0 lsds. A dividend rate strategy L = {ls} is called admissible if it satisfies that for a given M > 0, ls ∈ [0, M] is cadlag, adapted to the filtration {Fs} which is generated by {rs, Xs} and fulfilling Xs ≥ 0 up to the ruin time. We denote the set of admissible strategies by Uad[r, x]. Let L = {ls} be an admissible strategy and τ L denotes the ruin time of surplus process XL t with initial wealth X0 = x. The return function corresponding to L is V L (r, x) = E[ Z τ L 0 e −R s 0rudulsds], (r, x) ∈ R × {R+ ∪ 0}. (3.30) It means the dividend rate ls at time s is discounted by the factor e − R s 0 rudu. In the following, we write U y s as Us = R s 0  rudu with initial value r0 = y. Our target is to maximize the expected discounted dividends given the preference rate {rt}. We define value function as V (r, x) = sup L∈Uad[r,x] V L (r, x), (r, x) ∈ R × {R+ ∪ 0}. (3.31) The corresponding Hamilton-Jacobi-Bellman equation is [−(r + λ)V + a( ˆb − r)Vr + ˆδ 2 2 Vrr + cVx](r, x) + λ Z x 0 V (r, x − y)dG(y) + max 0≤l≤M l(1 − Vx(r, x)) = 0. (3.32) Given a continuously differentiable function ϕ(r, x) : R×[0, +∞) → R, we define the operator L[ϕ] :=[−(r + λ)ϕ + a( ˆb − r)ϕr + ˆδ 2 2 ϕrr + cϕx](r, x) + λ Z x 0 ϕ(r, x − y)dG(y) + max 0≤l≤M l(1 − ϕx(r, x)). (3.33) This definition will make us easier to state the definition and theorem. 3.1 Properties of the Value Function In this subsection we prove the boundedness and continuity of the value function V which is defined in (3.31). The continuity makes us easier to define viscosity solution. 

Lemma 3.1 The value function V is bounded.
Proof Via Fubini’s theorem, the value function satisfies
V (r, x) = sup L∈Uad[r,x] J L (r, x) ≤ E[ Z ∞ 0 e −Ur s M ds] = E[
Z ∞ 0 E[e −Ur s ]M ds]. (3.34)
Thanks to Borodin and Salminen (1998, p.525) [5], we can use the fact that E[e Ur s ] = e f(r,s) , where f(r, s) := −ˆbs + ˆδ 2 2σ 2 s − r − ˆb a (1 − e −as) + ˆδ 2 4a 3 1 − (2 − e −as) 2  . (3.35) 9



Sum them together and dividend by h, letting h → 0 and using the fact that l0 is arbitrarily, we obtain Lφ(r, x) ≤ 0. This proves the value function is viscosity supersolution of equation (3.32) .
Now we prove that value function is the viscosity subsolution of the corresponding HJB equation. Assume the contrary that there exists a point (r0, x0) that V is not the viscosity subsolution. By the definition of viscosity solution, there exists η > 0 and a continuously differentiable function ϕ 0 such that V (r0, x0) = ϕ 0 (r0, x0), ϕ 0 (r, x) ≥ V (r, x) for R×(0, +∞) and L[ϕ 0 ](r0, x0) = −2η < 0. Firstly, we assume that r0 ≥ 0 (r0 < 0 can be proved similarly). Consider function ϕˆ(r, x) = ϕ 0 (r, x) + η
x 2 0λ (x − x0) 2 + η λ (r − r0) 4 , (3.44) then we can notice that



ϕˆ(r0, x0) = ϕ(r0, x0);
ϕˆx(r0, x0) = ϕx(r0, x0);
ϕˆr(r0, x0) = ϕr(r0, x0);
ϕˆrr(r0, x0) = ϕrr(r0, x0).
And λ Z x0 0 ϕˆ(r0, x0 − y)βe−βydy =λ Z x0 0 [ϕ 0 (r0, x0 − y) + η x 2 0λ y 2 ]βe−βydy ≤λ Z x0 0 ϕ 0 (r0, x0 − y)βe−βydy + η.

We can get L[ ˆϕ](r0, x0) ≤ −η < 0.Since ˆϕ is nonnegative and continuously differentiable, we can find h ∈ (0,x2) such that
L[ ˆϕ](r, x) ≤ −η 2 < 0. (3.45) on (r, x) ∈ [r0−2h, r0+2h]×[x0−2h, x0+2h]. Let ψ be an nonnegative continuously differentiable
function with support included in (−1, 1) × (−1, 1) such that R 1
−1 R 1 −1 ψ(r, y)drdy = 1. We define νn : (−∞,∞) × [0, ∞) → R as the convolution νn(r, y) = 1 n2 Z Z√ |y−x| 2+|r−s| 2< 1 n ψ(n(r − s), n(y − x))(V (s, x) + ηh2 2λx2 0 + ηh4 2λ )dsdx. (3.46)
Since V is not defined on R−, in this integral we can extend V as V (r, y) = V (r, 0) + y for (r, y) ∈ R × (−∞, 0). We can find n0 large enough such that V (r, y) + ηh2 λx2 0 + ηh4 λ ≥ νn0 (r, y) ≥ V (r, y) + ηh2 4λx2 0 ηh4 4λ . (3.47) Let χ be a continuously differentiable function satisfying
13

## Concluding Remarks
In this paper we investigated the optimal dividend of insurance company under the assumption of stochastic interest rate and give out the explicit expression of the optimal strategy when interest rate follows geometric Brownian motion and claim size follows the exponential
distribution. For the case of Vasicek model, we didn’t give out the solution of value function but we explored the properties of the value function and we used the notion of viscosity solution to create the connection between value function and the HJB equation, which is important for the future study about the optimal strategy. 

When the discounting factor is given by geometric Brownian motion, we can see that the optimal strategy is still threshold strategy, except some changes in the parameters compared with the case of deterministic interest rate. This partly used the fact that surplus process is independent of the discounting factor, which provides a convenient condition for us to prove the optimality. We only considered the exponential claim in section 2. But we already started to explore more general case of claim distribution. We conjecture that in the setting of geometric Brownian motion, the optimal dividend is a band strategy if the claim follows a more general continuous distribution function G(y).

In section 3, we considered the dividend maximization problem when stochastic interest rate follows an Ornstein-Uhlenbeck Process. But we didn’t give out more regularities about value function. It’s quite hard to find an explicit expression of the dividend strategy. We will focus on optimal strategy in our future research. 

**Acknowledgements** This work is supported by the NSF of China (No. 11471171 and No. 11571189). 
Here I want to express my thanks to my supervisors Lihua Bai and Junyi Guo for their valuable insights and suggestions. 
I also want to say thanks to Jacques Rioux, Xiaoyi
Zhang. Thanks for their dedication to the improvement of this paper.

## References
[1] Albrecher, H., Thonhauser, S., Optimal dividend strategies for a risk process under force of interest, Insurance Math. Econom., 43 (2008), no. 1, 134-149.

[2] Albrecher, H., Thonhauser, S., Optimality results for dividend problems in insurance Rev. R. Acad. Cienc. Exactas Fs. Nat. Ser. A Math. RACSAM, 103 (2009), no. 2, 295-320.0

[3] Azcue, P., Muler, N., Optimal reinsurance and dividend distribution policies in the Cram´er-Lundberg model, Math. Finance, 15 (2005), no. 2, 261-308.

[4] Asmussen, S., Taksar, M., Controlled diffusion models for optimal dividend pay-out, Insurance Math. Econom., 20 (1997), no. 1, 1-15.

[5] Borodin, A. N., Salminen, P., Handbook of Brownian motion-facts and formulae, Birkh¨auser Verlag, Basel, 2002.

[6] Crandall, M. G., Ishii, H., User’s guide to viscosity solutions of second order partial differential equations, Bull. Amer. Math. Soc. (N.S.), 27 (1992), no. 1, 1-67.

[7] Crandall, M. G., Lions, P. L., Viscosity solutions of Hamilton-Jacobi equations, Trans. Amer. Math. Soc., 277 (1983), no. 1, 1-42.

[8] De Finetti B., Su un’impostazione alternativa della teoria collettiva del rischio, Transactions of the XVth international congress of Actuaries, Π 1957, no. 1, 433-443.

[9] Eisenberg, J., Optimal dividends under a stochastic interest rate, Insurance Math. Econom., 65 (2015), 259-266.

[10] Fleming, W. H., Soner, H. M., Controlled Markov processes and viscosity solutions, Second edition, Springer, New York, 2006.

[11] Gerber, H. U., Shiu, E. S. W., On optimal dividend strategies in the compound Poisson model, N. Am. Actuar. J., 10 (2006), no. 2, 76-93.

[12] Loeffen, R. L., On optimality of the barrier strategy in de Finetti’s dividend problem for spectrally negative Lvy processes, Ann. Appl. Probab., 18 (2008), no. 5, 1669-1680.

[13] Lions, P. L., Optimal control of diffusion processes and Hamilton-Jacobi-Bellman equations. I. The dynamic programming principle and applications, Comm. Partial Diff. Eqs., 8 (1983), no. 10, 1101-1174.

[14] Lions, P. L., Optimal control of diffusion processes and Hamilton-Jacobi-Bellman equations. II. Viscosity solutions and uniqueness., Comm. Partial Diff. Eqs., 8 (1983), no. 11,
1229-1276.

[15] Schmidli, H., Stochastic control in insurance, Springer, New York, 2008.

[16] Soner, H. M., Optimal control with state-space constraint. II., SIAM J. Control Optim., 24 (1986), no. 6, 1110-1122.

[17] Vasicek, O. A., An equilibrium characterization of the term structure, J. Financ. Econ., 5 (1977), no. 2, 177-188.

[18] Yong, J., Zhou, X. Y., Stochastic controls. Hamiltonian systems and HJB equations, Springer-Verlag, New York, 1999.

## Footnotes

[^1]: Not *"as a s geometric Brownian motion"* but **"as a geometric Brownian motion"**.
[^2]: Instead of *"up to the ruin time"* I would use *"up to ruin time"* or *"up to the time of ruin"*.
[^3]: Why **"explored"**? Why not **""explore"**?
[^4]: Instead of *"some properties of value function"* use *"some properties of **the** value function"*.
[^5]: Instead of *"we can not find the explicit expression of  the value function"*, use *"we can not find **an** explicit expression **for**  the value function"*
[^6]: *in the frame work* should be *in the* **framework**
[^7]: **Lundberge** should be **Lundberg**
[^8]: *"Later"* should be *"Later **,**"*
[^9]: Instead of *"influencing on the firm’s cost and profit"*, I would use either of  *"influencing **the** firm’s cost and profit"* or  *"with an influence on the firm’s cost and profit"*
[^10]: *"the interest **rates** is a function"* why is **rate** plural here?
[^11]: *"in the setting of surplus follows"* should be *"in the setting of surplus **following**"* or *"in the setting **where** surplus follows"*
[^12]: *"the claim size follows exponential distribution"* should be *"the claim size follows **an** exponential distribution"* 
[^13]: *"We explored dividend maximization problem"* should be *"We explored **the** dividend maximization problem"* 
[^14]:  *"We **explored** .... and find"* should be *"We **explore** .... and find"* or *"We **explored** .... and **found**"*
[^15]: Instead of *"but in the parameters"* I would use *"but the parameters do*".
[^16]: Instead of *"we consider Vasicek model, which the short rate"*, I would use *"we consider the Vasicek model, **for which** the short rate"*.
[^17]: Instead of **regularities** I think I would use **regularity properties**
[^18]: Not sure that  the word **"the"** is needed here.
[^19]:  *it do not need* should be *it **does** not need*
[^20]:  *merely requires the continuity of value function* should be *merely requires continuity of **the** value function*
