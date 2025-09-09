# Numerical study of the Bose–Hubbard model: 
# Exact Diagonalization & SDP  

This repository contains codes (in **Python** and **Julia**) to compute:  
- the **ground state energy** of the **Bose–Hubbard (BH) model**,
- the **one-body density matrix (OBDM)**,  
- the **two-body density matrix (TBDM)**  
- the **condensed fraction** of the ground state.  

---

## 📖 The Bose–Hubbard Model  

The Bose–Hubbard model originates from the **second quantization** of the Hamiltonian of interacting bosons in a periodic potential.  
In this model:  
- The periodic potential wells are assumed to be **deep enough** so that the space of possible positions for the bosons can be discretized into a **periodic lattice**.  
- Each **lattice site** corresponds to a potential well.  
- The Hamiltonian is then expressed in a **discrete Fock basis**, which encodes all possible distributions of bosons among the lattice sites.  

In a ring configuration of M sites, the Bose–Hubbard Hamiltonian reads:  

```math
\hat{H} = -J \sum_{i=1}^{M} (\hat{a}_i^\dagger \hat{a}_{i+1} + \hat{a}_{i+1}^\dagger \hat{a}_{i}) 
+ \frac{U}{2} \sum_{i=1}^{M} \hat{n}_i (\hat{n}_i - 1) 
```
where  

$$J$$ : hopping amplitude  

$$U$$ : on-site interaction strength 

$$\hat{a}_i^\dagger, \hat{a}_i$$ : bosonic creation/annihilation operators  

$$\hat{n}_i = \hat{a}_i^\dagger \hat{a}_i$$ : number operator


## ⚙️ Implemented Methods  

These codes implement and compare two approaches to compute the **ground state properties** (energy, OBDM, TBDM):  

1. **Exact Diagonalization (ED)**  
   - Construction of the Hamiltonian in the **Fock basis**.  
   - Direct diagonalization of the Hamiltonian matrix.  

2. **Semidefinite Programming (SDP) approach**  
   - Uses the **OBDM** and **TBDM** as optimization variables.  
   - Solved using the **CSD solver** in Julia.  
   - Allows access to larger system sizes beyond the reach of ED.  

---

## Python Code  

- Contains **only the classical ED method**.  
- Improved efficiency thanks to the **Ponomarev** ordering function, which optimizes the indexing of Fock states compared to the naive construction.  

---

## Julia Code  

- Contains both:  
  - **Exact Diagonalization**,  
  - **SDP approach** (via the **CSD solver**).  
- Provides a more scalable framework for studying larger Bose–Hubbard systems.  

---

## For more infos: 

These codes were developed during a **research internship** supervised by **Mr. Bruno Julia** at **Barcelona University**.  
A **poster presentation** of the internship is included in this repository. 
