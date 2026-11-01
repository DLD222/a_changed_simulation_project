# SPDX-License-Identifier: BSD-3-Clause
import argparse, math, os, re, subprocess, sys
import numpy as np
C=299792458.0
BIN="build_2d/bin/warpx.2d.MPI.OMP.DP.PDP.EB"
TMPL="""max_step = 900
amr.n_cell = {n} {n}
amr.max_grid_size = {n}
amr.max_level = 0
geometry.dims = 2
geometry.prob_lo = -2.0 -2.0
geometry.prob_hi =  2.0  2.0
warpx.cfl = 0.4
algo.maxwell_solver = ect
algo.ect_conformal_E = 2
algo.particle_shape = 1
algo.em_solver_medium = macroscopic
macroscopic.epsilon = {epsr}*epsilon0
macroscopic.mu = mu0
macroscopic.sigma = 0
boundary.field_lo = pec pec
boundary.field_hi = pec pec
warpx.eb_implicit_function = "z - 100"
warpx.E_ext_grid_init_style = parse_E_ext_grid_function
warpx.Ex_external_grid_function(x,y,z) = "0"
warpx.Ey_external_grid_function(x,y,z) = "exp(-((z+1.0)/{w})^2)*cos({k0}*(z-({zc})))"
warpx.Ez_external_grid_function(x,y,z) = "0"
warpx.B_ext_grid_init_style = parse_B_ext_grid_function
warpx.Bx_external_grid_function(x,y,z) = "-exp(-((z+1.0)/{w})^2)*cos({k0}*(z-({zc})))/{v}"
warpx.By_external_grid_function(x,y,z) = "0"
warpx.Bz_external_grid_function(x,y,z) = "0"
warpx.reduced_diags_names = p1 p2 p3 pE
p1.type = FieldProbe
p1.intervals = 1
p1.probe_geometry = Point
p1.x_probe = 0.0
p1.z_probe = -0.5
p2.type = FieldProbe
p2.intervals = 1
p2.probe_geometry = Point
p2.x_probe = 0.0
p2.z_probe = 0.0
p3.type = FieldProbe
p3.intervals = 1
p3.probe_geometry = Point
p3.x_probe = 0.0
p3.z_probe = 0.5
pE.type = FieldEnergy
pE.intervals = 1
"""
def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("run"); ap.add_argument("--epsr",type=float,default=2.25)
    ap.add_argument("--n",type=int,default=128); ap.add_argument("--w",type=float,default=0.5)
    ap.add_argument("--lam",type=float,default=1.0); ap.add_argument("--zc",type=float,default=-1.0)
    ap.add_argument("--make",action="store_true"); ap.add_argument("--execute",action="store_true")
    a=ap.parse_args()
    k0=2*math.pi/a.lam; v=C/math.sqrt(a.epsr)
    if a.make:
        os.makedirs(a.run,exist_ok=True)
        open(os.path.join(a.run,"inputs"),"w").write(TMPL.format(n=a.n,epsr=a.epsr,w=a.w,k0=k0,zc=a.zc,v=v))
    if a.execute:
        env=dict(os.environ,OMP_NUM_THREADS="8")
        with open(os.path.join(a.run,"run.log"),"w") as f:
            subprocess.run([BIN,"inputs"],cwd=a.run,stdout=f,stderr=subprocess.STDOUT,env=env,check=False)
    # DFT phase speed
    dx=4.0/a.n; dt=0.4*dx/(C*math.sqrt(2)); omega=k0*C
    zs=[]; ph=[]
    for p in ("p1","p2","p3"):
        d=np.loadtxt(os.path.join(a.run,"diags/reducedfiles",p+".txt"),skiprows=1)
        t=d[:,1]; ey=d[:,6]
        A=np.sum(ey*np.exp(1j*omega*t))
        zs.append(float(d[0,4])); ph.append(np.angle(A))
    zs=np.array(zs); ph=np.unwrap(np.array(ph))
    k=np.polyfit(zs,ph,1)[0]; v_num=omega/abs(k)
    # material-weighted energy
    e=np.loadtxt(os.path.join(a.run,"diags/reducedfiles/pE.txt"),skiprows=1)
    emat=a.epsr*e[:,3]+e[:,4]
    # old arrival-time estimate at p1
    d=np.loadtxt(os.path.join(a.run,"diags/reducedfiles/p1.txt"),skiprows=1)
    ey=np.abs(d[:,6]); m=ey.max(); st=int(d[int(np.argmax(ey>0.05*m)),0])
    v_arr=0.5/(st*dt) if st>0 else float("nan")
    print(f"epsr={a.epsr:5.2f}  v_theory={v/C:.4f}c  v_dft={v_num/C:.4f}c  rel={v_num/v-1:+.2%}  "
          f"v_arrival={v_arr/C:.4f}c  energy(material) final/init={emat[-1]/emat[0]:.5f} max/min={emat.max()/emat.min():.4f}")
if __name__=="__main__": main()
