"""Movimento não relativístico em campos elétrico e magnético uniformes (SI)."""

import argparse
import json
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp


@dataclass(frozen=True)
class Parametros:
    q_m: float = 1.0  # C/kg; mantém o sinal da carga
    E: tuple = (0.0, 0.2, 0.0)  # V/m
    B: tuple = (0.0, 0.0, 1.0)  # T
    r0: tuple = (0.0, 0.0, 0.0)  # m
    v0: tuple = (1.0, 0.0, 0.2)  # m/s
    duracao: float = 30.0  # s
    pontos: int = 600

    def validar(self):
        for nome in ("E", "B", "r0", "v0"):
            vetor = np.asarray(getattr(self, nome), dtype=float)
            if vetor.shape != (3,) or not np.all(np.isfinite(vetor)):
                raise ValueError(f"{nome} deve conter três números finitos.")
        if not np.isfinite(self.q_m):
            raise ValueError("q/m deve ser finito.")
        if not np.isfinite(self.duracao) or self.duracao <= 0:
            raise ValueError("A duração deve ser positiva e finita.")
        if not isinstance(self.pontos, (int, np.integer)) or self.pontos < 2:
            raise ValueError("Use pelo menos dois pontos.")


def derivadas(t, estado, q_m, E, B):
    """Retorna [vx, vy, vz, ax, ay, az], nunca a força isolada."""
    vx, vy, vz = estado[3:]
    Ex, Ey, Ez = E
    Bx, By, Bz = B
    return [vx, vy, vz,
            q_m * (Ex + vy * Bz - vz * By),
            q_m * (Ey + vz * Bx - vx * Bz),
            q_m * (Ez + vx * By - vy * Bx)]


def solucao_exata(p, tempos):
    """Referência analítica para campos constantes, inclusive B=0 e q/m=0."""
    t = np.asarray(tempos, dtype=float)[None, :]
    E, B, r0, v0 = (np.asarray(a, dtype=float) for a in (p.E, p.B, p.r0, p.v0))
    modulo_B = np.linalg.norm(B)
    if modulo_B == 0 or p.q_m == 0:
        a = p.q_m * E[:, None]
        r = r0[:, None] + v0[:, None] * t + 0.5 * a * t**2
        v = v0[:, None] + a * t
        return np.vstack((r, v))

    b = B / modulo_B
    omega = p.q_m * modulo_B  # frequência com sinal
    deriva = np.cross(E, B) / modulo_B**2
    E_paralelo = np.dot(E, b) * b
    u = v0 - deriva
    u_paralelo = np.dot(u, b) * b
    u_perpendicular = u - u_paralelo
    rotacao = np.cross(u_perpendicular, b)
    angulo = omega * t
    # sinc evita perda de precisão em sin(omega*t)/omega quando omega*t é pequeno.
    integral_seno = t * np.sinc(angulo / np.pi)
    integral_cosseno = 0.5 * omega * t**2 * np.sinc(angulo / (2 * np.pi))**2
    a_paralela = p.q_m * E_paralelo[:, None]
    r = (r0[:, None] + (deriva + u_paralelo)[:, None] * t
         + 0.5 * a_paralela * t**2
         + u_perpendicular[:, None] * integral_seno
         + rotacao[:, None] * integral_cosseno)
    v = (deriva[:, None] + u_paralelo[:, None] + a_paralela * t
         + u_perpendicular[:, None] * np.cos(angulo)
         + rotacao[:, None] * np.sin(angulo))
    return np.vstack((r, v))


def simular(p, metodo="RK45", rtol=1e-8, atol=1e-10):
    p.validar()
    omega = abs(p.q_m) * np.linalg.norm(p.B)
    periodo = 2 * np.pi / omega if omega > 0 else np.inf
    # Pelo menos 80 amostras por volta para desenhar a órbita sem serrilhamento.
    quantidade = max(p.pontos, int(np.ceil(80 * p.duracao / periodo)) + 1)
    if quantidade > 200_000:
        raise ValueError("Simulação muito extensa: reduza a duração ou a frequência ciclotrônica.")
    tempos = np.linspace(0, p.duracao, quantidade)
    sol = solve_ivp(
        derivadas, (0, p.duracao), [*p.r0, *p.v0],
        args=(p.q_m, p.E, p.B), method=metodo, t_eval=tempos,
        rtol=rtol, atol=atol, max_step=periodo / 20,
    )
    if not sol.success:
        raise RuntimeError(sol.message)
    return sol


def diagnosticos(p, sol):
    referencia = solucao_exata(p, sol.t)
    r, v = sol.y[:3], sol.y[3:]
    cinetica = 0.5 * np.sum(v**2, axis=0)  # J/kg
    potencial = -p.q_m * np.asarray(p.E) @ r  # J/kg; zero na origem
    total = cinetica + potencial
    erro = np.linalg.norm(r - referencia[:3], axis=0)
    return referencia, cinetica, potencial, total, erro


def criar_painel(p, sol):
    import matplotlib.pyplot as plt

    referencia, K, U, H, erro = diagnosticos(p, sol)
    r = sol.y[:3]
    fig = plt.figure(figsize=(15, 9), layout="constrained")
    fig.suptitle("Força de Lorentz · simulação e verificação física", fontsize=17)
    ax = fig.add_subplot(231, projection="3d")
    ax.plot(*r, label="Numérica", color="crimson")
    ax.plot(*referencia[:3, ::10], ".", ms=2, label="Exata", color="black")
    ax.scatter(*r[:, 0], color="green", label="Início")
    ax.set(xlabel="X (m)", ylabel="Y (m)", zlabel="Z (m)", title="Trajetória 3D")
    ajustar_escala_3d(ax, r)
    ax.legend(fontsize=8)

    ax = fig.add_subplot(232)
    ax.plot(r[0], r[1], color="crimson")
    ax.set(xlabel="X (m)", ylabel="Y (m)", title="Projeção XY")
    ax.set_aspect("equal", adjustable="datalim")

    ax = fig.add_subplot(233)
    for valores, nome in zip(r, "xyz"):
        ax.plot(sol.t, valores, label=nome)
    ax.set(xlabel="Tempo (s)", ylabel="Posição (m)", title="Posição ao longo do tempo")
    ax.legend()

    ax = fig.add_subplot(234)
    ax.plot(sol.t, np.linalg.norm(sol.y[3:], axis=0))
    ax.set(xlabel="Tempo (s)", ylabel="Rapidez (m/s)", title="Módulo da velocidade")

    ax = fig.add_subplot(235)
    for valores, nome in ((K, "Cinética"), (U, "Potencial elétrica"), (H, "Total")):
        ax.plot(sol.t, valores, label=nome)
    ax.set(xlabel="Tempo (s)", ylabel="Energia / massa (J/kg)", title="Conservação da energia")
    ax.legend(fontsize=8)

    ax = fig.add_subplot(236)
    ax.plot(sol.t, erro)
    ax.set(xlabel="Tempo (s)", ylabel="Erro de posição (m)", title="Distância à solução exata")
    for ax in fig.axes[1:]:
        ax.grid(alpha=0.25)
    return fig


def ajustar_escala_3d(ax, r):
    amplitudes = np.ptp(r, axis=1)
    margem = max(float(amplitudes.max()) * 0.06, 0.05)
    limites = np.array([(a.min() - margem, a.max() + margem) for a in r])
    ax.set_xlim(*limites[0])
    ax.set_ylim(*limites[1])
    ax.set_zlim(*limites[2])
    ax.set_box_aspect(limites[:, 1] - limites[:, 0])


def criar_animacao(p, sol, fps=30, max_quadros=360):
    import matplotlib.pyplot as plt
    from matplotlib.animation import FuncAnimation

    fig = plt.figure(figsize=(9, 7), layout="constrained")
    ax = fig.add_subplot(111, projection="3d")
    r = sol.y[:3]
    ajustar_escala_3d(ax, r)
    ax.set(xlabel="X (m)", ylabel="Y (m)", zlabel="Z (m)", title="Partícula sob a força de Lorentz")
    ax.plot(*r, color="gray", alpha=0.2, label="Percurso completo")
    linha, = ax.plot([], [], [], color="crimson", label="Trajetória")
    ponto, = ax.plot([], [], [], "o", color="royalblue", label="Partícula")
    texto = ax.text2D(0.03, 0.95, "", transform=ax.transAxes)
    ax.legend(loc="lower left")
    indices = np.unique(np.linspace(0, len(sol.t) - 1, min(max_quadros, len(sol.t))).astype(int))

    def atualizar(i):
        linha.set_data_3d(*r[:, :i + 1])
        ponto.set_data_3d(*r[:, i:i + 1])
        texto.set_text(f"Tempo físico: {sol.t[i]:.2f} s")
        return linha, ponto, texto

    animacao = FuncAnimation(fig, atualizar, frames=indices, interval=1000 / fps,
                             blit=False, repeat=True, cache_frame_data=False)
    return fig, animacao


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cenario", choices=["cruzados", "livre", "eletrico", "magnetico"], default="cruzados")
    parser.add_argument("--q-m", type=float, default=1.0, help="Razão carga/massa em C/kg")
    parser.add_argument("--E", nargs=3, type=float, metavar=("EX", "EY", "EZ"))
    parser.add_argument("--B", nargs=3, type=float, metavar=("BX", "BY", "BZ"))
    parser.add_argument("--v0", nargs=3, type=float, default=(1, 0, 0.2))
    parser.add_argument("--r0", nargs=3, type=float, default=(0, 0, 0))
    parser.add_argument("--duracao", type=float, default=30.0)
    parser.add_argument("--pontos", type=int, default=600, help="Número mínimo de amostras")
    parser.add_argument("--metodo", choices=["RK45", "DOP853"], default="RK45")
    parser.add_argument("--saida", type=Path, default=Path("output"))
    parser.add_argument("--no-show", action="store_true", help="Gera arquivos sem abrir janelas")
    parser.add_argument("--gif", action="store_true", help="Também exporta a animação GIF")
    args = parser.parse_args()
    campos = {
        "cruzados": ((0, 0.2, 0), (0, 0, 1)),
        "livre": ((0, 0, 0), (0, 0, 0)),
        "eletrico": ((0, 0.2, 0), (0, 0, 0)),
        "magnetico": ((0, 0, 0), (0, 0, 1)),
    }
    E, B = campos[args.cenario]
    p = Parametros(q_m=args.q_m, E=tuple(args.E) if args.E is not None else E,
                   B=tuple(args.B) if args.B is not None else B,
                   r0=tuple(args.r0), v0=tuple(args.v0),
                   duracao=args.duracao, pontos=args.pontos)
    try:
        sol = simular(p, args.metodo)
    except (ValueError, RuntimeError) as exc:
        parser.error(str(exc))

    import matplotlib
    if args.no_show:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    args.saida.mkdir(parents=True, exist_ok=True)
    ref, K, U, H, erro = diagnosticos(p, sol)
    painel = criar_painel(p, sol)
    painel.savefig(args.saida / "diagnosticos.png", dpi=150)
    dados = np.column_stack((sol.t, sol.y.T, K, U, H, erro))
    np.savetxt(args.saida / "trajetoria.csv", dados, delimiter=",", comments="",
               header="t_s,x_m,y_m,z_m,vx_m_s,vy_m_s,vz_m_s,K_J_kg,U_J_kg,H_J_kg,erro_posicao_m")
    omega = abs(p.q_m) * np.linalg.norm(p.B)
    resumo = {
        "parametros_SI": {**p.__dict__, "pontos_efetivos": len(sol.t)},
        "metodo": args.metodo,
        "avaliacoes_da_equacao": sol.nfev,
        "erro_maximo_posicao_m": float(erro.max()),
        "variacao_maxima_energia_especifica_J_kg": float(np.max(np.abs(H - H[0]))),
        "periodo_ciclotronico_s": float(2 * np.pi / omega) if omega else None,
        "deriva_E_cross_B_m_s": (np.cross(p.E, p.B) / np.dot(p.B, p.B)).tolist()
        if np.linalg.norm(p.B) > 0 and p.q_m != 0 else None,
    }
    (args.saida / "resumo.json").write_text(json.dumps(resumo, indent=2, ensure_ascii=False) + "\n")
    print(json.dumps(resumo, indent=2, ensure_ascii=False))
    print(f"Arquivos gravados em: {args.saida.resolve()}")
    animacao = None  # conserva a referência enquanto a janela estiver aberta
    if args.gif or not args.no_show:
        figura_animada, animacao = criar_animacao(p, sol)
        if args.gif:
            animacao.save(args.saida / "animacao.gif", writer="pillow", fps=30, dpi=85)
    if not args.no_show:
        plt.show()
    plt.close("all")


if __name__ == "__main__":
    main()
