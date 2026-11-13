"""Verificações físicas: execute python -m unittest -v."""

import unittest
from dataclasses import replace

import numpy as np

from simulacao import Parametros, derivadas, diagnosticos, simular, solucao_exata


class TestFisica(unittest.TestCase):
    def test_movimento_livre(self):
        p = Parametros(E=(0, 0, 0), B=(0, 0, 0))
        sol = simular(p)
        esperado = np.asarray(p.r0)[:, None] + np.asarray(p.v0)[:, None] * sol.t
        np.testing.assert_allclose(sol.y[:3], esperado, atol=1e-11)

    def test_campo_eletrico_aceleracao_constante(self):
        p = Parametros(B=(0, 0, 0))
        sol = simular(p)
        np.testing.assert_allclose(sol.y[1], 0.1 * sol.t**2, atol=1e-10)
        np.testing.assert_allclose(sol.y[4], 0.2 * sol.t, atol=1e-10)

    def test_campo_magnetico_preserva_rapidez(self):
        p = Parametros(E=(0, 0, 0))
        sol = simular(p)
        rapidez2 = np.sum(sol.y[3:]**2, axis=0)
        np.testing.assert_allclose(rapidez2, np.dot(p.v0, p.v0), atol=2e-7)

    def test_cruzados_concordam_com_formula_do_relatorio(self):
        p = Parametros()
        sol = simular(p)
        t = sol.t
        esperado = np.array([0.2*t + 0.8*np.sin(t), 0.8*(np.cos(t)-1), 0.2*t])
        np.testing.assert_allclose(sol.y[:3], esperado, atol=2e-7)
        _, _, _, H, _ = diagnosticos(p, sol)
        np.testing.assert_allclose(H, H[0], atol=2e-7)

    def test_carga_negativa_inverte_forca(self):
        estado = [0, 0, 0, 1, 0, 0]
        positiva = derivadas(0, estado, 1, (0, 0, 0), (0, 0, 1))
        negativa = derivadas(0, estado, -1, (0, 0, 0), (0, 0, 1))
        self.assertEqual(positiva[4], -1)
        self.assertEqual(negativa[4], 1)

    def test_campos_obliquos_e_sinal_da_carga(self):
        for q_m in (1.3, -0.7):
            with self.subTest(q_m=q_m):
                p = Parametros(q_m=q_m, E=(0.3, -0.1, 0.5), B=(0.2, 0.5, 0.8),
                               r0=(1, 2, -1), v0=(0.2, -0.4, 0.8))
                sol = simular(p)
                np.testing.assert_allclose(sol.y, solucao_exata(p, sol.t), atol=2e-7)

    def test_particula_neutra(self):
        p = Parametros(q_m=0)
        sol = simular(p)
        np.testing.assert_allclose(sol.y, solucao_exata(p, sol.t), atol=1e-11)

    def test_periodo_retorna_ao_mesmo_ponto_no_plano(self):
        p = Parametros(E=(0, 0, 0), duracao=2*np.pi)
        sol = simular(p)
        np.testing.assert_allclose(sol.y[:2, -1], p.r0[:2], atol=1e-7)

    def test_parametros_invalidos(self):
        for p in (replace(Parametros(), duracao=-1), replace(Parametros(), pontos=1),
                  replace(Parametros(), E=(0, float("nan"), 0))):
            with self.assertRaises(ValueError):
                simular(p)


if __name__ == "__main__":
    unittest.main()
