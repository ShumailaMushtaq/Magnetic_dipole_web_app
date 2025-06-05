import streamlit as st
import numpy as np
import matplotlib.pyplot as plt
import magpylib as magpy
from magpylib import misc
import os           

st.set_page_config(layout="centered")
st.title("Magnetic Dipole Field")

st.markdown(
    """
    **Instructions:**  
    1. Enter your **Peak Drive Voltage** \(V_0\), coil **Resistance** \(R\), coil **Inductance** \(L\),  
       coil **Diameter** \(d\), **Number of Turns** \(N\), and a small **axial offset** \(z_0\) in the sidebar.  
    2. Click **Compute & Plot**.  
    3. You’ll see:  
       - The **sinusoidal coil current**, denoted $I_{\mathrm{coil}}(t)$, i.e.\  
         $$I_{\mathrm{coil}}(t) = I_0\,\sin(\omega t)\,,\quad \omega = 2\pi f.$$
       - The magnetic‐dipole field magnitude $|\mathbf{B}(t)|$ at \((z_0,\,0,\,0)\).  
       - The computed **peak current** $I_0$, which is the amplitude of that sine wave.  
    4. When you’re done, click **Exit App** to stop the server without needing Ctrl+C.
    """
)

# ── Sidebar Inputs ─────────────────────────────────────────────────────────────
st.sidebar.header("Coil & Drive Parameters")

V0 = st.sidebar.number_input(
    "Peak Drive Voltage $V_0$ (V)",
    value=10.0,
    step=1.0,
    format="%.3f"
)
f = st.sidebar.number_input(
    "Frequency $f$ (Hz)",
    value=50.0,
    step=1.0,
    format="%.3f"
)
R_coil = st.sidebar.number_input(
    "Coil Resistance $R$ (Ω)",
    value=1.0,
    step=0.1,
    format="%.3f"
)
L_coil = st.sidebar.number_input(
    "Coil Inductance $L$ (H)",
    value=1e-3,
    step=1e-4,
    format="%.6f"
)
diameter = st.sidebar.number_input(
    "Coil Diameter $d$ (m)",
    value=0.10,
    step=0.01,
    format="%.3f"
)
num_turns = st.sidebar.number_input(
    "Number of Turns $N$",
    value=1,
    step=1,
    min_value=1
)
z0 = st.sidebar.number_input(
    "Axial Offset $z_0$ (m)",
    value=0.01,
    step=0.005,
    format="%.3f"
)
samples_per_period = st.sidebar.number_input(
    "Samples per Period",
    value=500,
    step=100,
    min_value=100
)

compute_button = st.sidebar.button("Compute & Plot")
exit_button = st.sidebar.button("Exit App")  # <— new

# ── Exit Logic ─────────────────────────────────────────────────────────────────
if exit_button:
    st.write("Shutting down the app…")
    os._exit(0)  # immediately kill the Python process

# ── Main Computation & Plot ────────────────────────────────────────────────────
if compute_button:
    # 1. Compute peak current amplitude I0 = V0 / |Z|, where |Z| = sqrt(R^2 + (ωL)^2)
    omega = 2 * np.pi * f
    Z_mag = np.sqrt(R_coil**2 + (omega * L_coil)**2)
    I0 = V0 / Z_mag

    # 2. Time array over one period T = 1/f
    T = 1.0 / f
    t = np.linspace(0, T, int(samples_per_period))

    # 3. Sinusoidal coil current: I_coil(t) = I0 * sin(ω t)
    I_coil = I0 * np.sin(omega * t)

    # 4. Coil area A = π (d/2)^2
    radius = diameter / 2.0
    A_tx = np.pi * radius**2

    # 5. Create a “dummy” dipole—moment overwritten each loop iteration
    dipole = misc.Dipole(
        moment=(0.0, 0.0, num_turns * A_tx),  # placeholder
        position=(0.0, 0.0, 0.0),
    )

    # 6. Loop: set dipole moment = I_coil(t)*N*A, compute B at (z0, 0, 0)
    B_vals = np.empty_like(t)
    for i, I_inst in enumerate(I_coil):
        m_val = I_inst * num_turns * A_tx
        dipole.moment = (0.0, 0.0, m_val)
        B_vec = magpy.getB(dipole, (z0, 0.0, 0.0))
        B_vals[i] = np.linalg.norm(B_vec)

    # 7. Plot I_coil(t) and |B|(t)
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(6, 5), tight_layout=True)

    ax1.plot(t, I_coil, color="C0")
    ax1.set_title("Instantaneous Coil Current $I_{\\mathrm{coil}}(t)$")
    ax1.set_xlabel("Time (s)")
    ax1.set_ylabel("$I_{\\mathrm{coil}}(t)$ [A]")
    ax1.grid(True)

    ax2.plot(t, B_vals, color="C1")
    ax2.set_title(r"Magnetic Dipole Field $|\mathbf{B}(t)|$ at $(z_0,0,0)$")
    ax2.set_xlabel("Time (s)")
    ax2.set_ylabel(r"$|\mathbf{B}(t)|$ [T]")
    ax2.grid(True)

    st.pyplot(fig)

    # 8. Display computed peak current I0
    st.markdown(f"**Computed peak current $I_0$:** {I0:.4f} A")

    st.stop()  # <— prevents any further execution once results are shown
