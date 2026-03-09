import numpy as np
import matplotlib.pyplot as plt


def matrice_rotation(theta_1, theta_2, theta_3, theta_4, theta_5, theta_6):
    wRa = np.array([
        [np.cos(theta_1), 0, np.sin(theta_1)],
        [0, 1, 0],
        [-np.sin(theta_1), 0, np.cos(theta_1)]
    ])

    aRb = np.array([
        [np.cos(theta_2), -np.sin(theta_2), 0],
        [np.sin(theta_2),  np.cos(theta_2), 0],
        [0, 0, 1]
    ])

    bRc = np.array([
        [np.cos(theta_3), -np.sin(theta_3), 0],
        [np.sin(theta_3),  np.cos(theta_3), 0],
        [0, 0, 1]
    ])

    cRd = np.array([
        [1, 0, 0],
        [0,  np.cos(theta_4), -np.sin(theta_4)],
        [0,  np.sin(theta_4),  np.cos(theta_4)]
    ])

    dRe = np.array([
        [np.cos(theta_5), -np.sin(theta_5), 0],
        [np.sin(theta_5),  np.cos(theta_5), 0],
        [0, 0, 1]
    ])

    eRt = np.array([
        [1, 0, 0],
        [0,  np.cos(theta_6), -np.sin(theta_6)],
        [0,  np.sin(theta_6),  np.cos(theta_6)]
    ])

    return wRa, aRb, bRc, cRd, dRe, eRt


def analyse_tranche():
    # Angles (rad)
    theta_1, theta_2, theta_3, theta_4, theta_5, theta_6 = 0, -0.3, 0, 0, 0.5, -1.6

    # Dimensions (mètres)
    hg, hd, lb = 0.1, 0.05, 0.15
    lx_cam, ly_cam = 0.8, 0.7
    lpx, lpz = 0.5994, 0.1991

    r_w_VwW = np.array([[lx_cam], [ly_cam], [0]])

    # Calcul de r_p_TwP avec configuration de la figure 2
    _, r_w_TwW_figure2 = cinematique(-0.4, -1.2, 0, 0, -0.3708, 0)
    r_w_PwW = np.array([[lpx], [0], [lpz]])
    r_w_PwT = r_w_PwW - r_w_TwW_figure2
    pRw_figure2 = calcul_pRw(-0.4, -1.2, 0, 0, -0.3708, 0)
    r_p_PwT = pRw_figure2 @ r_w_PwT
    # print(f"r_p_TwP: {r_p_PwT}")

    # Coord TIi dans le repère v
    r_v_TI1wV = np.array([ [0.158920], [0.013914], [0.028686] ])
    r_v_TI2wV = np.array([ [0.157470], [0.021067], [0.008891] ])
    r_v_TI3wV = np.array([ [0.153781], [0.039266], [-0.040587] ])
    r_v_TI4wV = np.array([ [0.152420], [0.04597], [-0.060395] ])
    r_v_TI5wV = np.array([ [0.150931], [0.053326], [-0.080185] ])

    # Coord DFi dans le repère v
    r_v_DF1wV = np.array([ [0.153758], [0.039379], [-0.020575] ])
    r_v_DF2wV = np.array([ [0.145698], [0.079138], [-0.039398] ])
    r_v_DF3wV = np.array([ [0.153932], [0.038521], [0.009411] ])
    r_v_DF4wV = np.array([ [0.152097], [0.047573], [0.035692] ])
    r_v_DF5wV = np.array([ [0.146104], [0.077134], [0.030571] ])

    # Matrice de rotation
    wRv = np.array([
        [-1, 0, 0],
        [0, -1, 0],
        [0,  0, 1]
    ])

    pRt = np.array([
        [0, 1, 0],
        [0, 0, 1],
        [1, 0, 0]
    ])

    wRa, aRb, bRc, cRd, dRe, eRt = matrice_rotation(
        theta_1, theta_2, theta_3, theta_4, theta_5, theta_6
    )
    aRw = wRa.T
    bRa = aRb.T
    cRb = bRc.T
    dRc = cRd.T
    eRd = dRe.T
    tRe = eRt.T

    pRv = pRt @ tRe @ eRd @ dRc @ cRb @ bRa @ aRw @ wRv
    pRw = pRt @ tRe @ eRd @ dRc @ cRb @ bRa @ aRw

    # TI, DF dans le repère p
    r_p_TI1wV = pRv @ r_v_TI1wV
    r_p_TI2wV = pRv @ r_v_TI2wV
    r_p_TI3wV = pRv @ r_v_TI3wV
    r_p_TI4wV = pRv @ r_v_TI4wV
    r_p_TI5wV = pRv @ r_v_TI5wV

    r_p_DF1wV = pRv @ r_v_DF1wV
    r_p_DF2wV = pRv @ r_v_DF2wV
    r_p_DF3wV = pRv @ r_v_DF3wV
    r_p_DF4wV = pRv @ r_v_DF4wV
    r_p_DF5wV = pRv @ r_v_DF5wV

    # Transformation du vecteur caméra
    r_p_VwW = pRw @ r_w_VwW

    # Calcul de la position w -> t
    _, r_w_TwW = cinematique(theta_1, theta_2, theta_3, theta_4, theta_5, theta_6)
    r_p_TwW = pRw @ r_w_TwW

    # TIi, repère p || p -> TIi
    r_p_TI1wP = (r_p_TI1wV + r_p_VwW) - (r_p_PwT + r_p_TwW)
    r_p_TI2wP = (r_p_TI2wV + r_p_VwW) - (r_p_PwT + r_p_TwW)
    r_p_TI3wP = (r_p_TI3wV + r_p_VwW) - (r_p_PwT + r_p_TwW)
    r_p_TI4wP = (r_p_TI4wV + r_p_VwW) - (r_p_PwT + r_p_TwW)
    r_p_TI5wP = (r_p_TI5wV + r_p_VwW) - (r_p_PwT + r_p_TwW)

    # DFi, repère p || p -> DFi
    r_p_DF1wP = (r_p_DF1wV + r_p_VwW) - (r_p_PwT + r_p_TwW)
    r_p_DF2wP = (r_p_DF2wV + r_p_VwW) - (r_p_PwT + r_p_TwW)
    r_p_DF3wP = (r_p_DF3wV + r_p_VwW) - (r_p_PwT + r_p_TwW)
    r_p_DF4wP = (r_p_DF4wV + r_p_VwW) - (r_p_PwT + r_p_TwW)
    r_p_DF5wP = (r_p_DF5wV + r_p_VwW) - (r_p_PwT + r_p_TwW)

    # angle theo
    delta_h = hg - hd
    theta_theorique = -np.degrees(np.arctan(delta_h / lb))      # voir fig 1 pour -

    ti_x = np.array([r_p_TI1wP[0,0], r_p_TI2wP[0,0], r_p_TI3wP[0,0], r_p_TI4wP[0,0], r_p_TI5wP[0,0]])
    ti_y = np.array([r_p_TI1wP[1,0], r_p_TI2wP[1,0], r_p_TI3wP[1,0], r_p_TI4wP[1,0], r_p_TI5wP[1,0]])

    # y = ax + b
    A = np.column_stack([ti_x, np.ones(len(ti_x))])
    a_b, _, _, _ = np.linalg.lstsq(A, ti_y, rcond=None)
    a, b = a_b[0], a_b[1]

    theta_measure = np.degrees(np.arctan(a))

    # Vérification 0.2°
    if np.abs(theta_theorique - theta_measure) <= 0.2:
        print(f"différence d'angle: {theta_theorique - theta_measure}\nTest de 0.2° réussi")
    else:
        print(f"différence d'angle: {theta_theorique - theta_measure}\nTest de 0.2° échoué")

    # Vérification ±1 mm
    # print(b)    # 0.10125892090668703
    if np.abs(b - hg) <= 0.001:
        print(f"différence d'hauteur: {b - hg}\nTest de ±1 mm réussi")
    else:
        print(f"différence d'hauteur: {b - hg}\nTest de ±1 mm échoué")

    # DF
    matrice_r_p_DFiwP = np.stack((r_p_DF1wP, r_p_DF2wP, r_p_DF3wP, r_p_DF4wP, r_p_DF5wP))

    for i in range(len(matrice_r_p_DFiwP)):
        x = matrice_r_p_DFiwP[i, 0]
        y = matrice_r_p_DFiwP[i, 1]

        if verifie_zone(x, y):
            print(f"Défaut détecté: DF{i}: ({x}, {y})")

    df_x = np.array([r_p_DF1wP[0,0], r_p_DF2wP[0,0], r_p_DF3wP[0,0], r_p_DF4wP[0,0], r_p_DF5wP[0,0]])
    df_y = np.array([r_p_DF1wP[1,0], r_p_DF2wP[1,0], r_p_DF3wP[1,0], r_p_DF4wP[1,0], r_p_DF5wP[1,0]])

    dessiner_zone(lb, hg, hd, df_x, df_y, ti_x, ti_y)


def calcul_pRw(theta_1, theta_2, theta_3, theta_4, theta_5, theta_6):
    wRa, aRb, bRc, cRd, dRe, eRt = matrice_rotation(
        theta_1, theta_2, theta_3, theta_4, theta_5, theta_6
    )
    aRw = wRa.T
    bRa = aRb.T
    cRb = bRc.T
    dRc = cRd.T
    eRd = dRe.T
    tRe = eRt.T

    pRt = np.array([
        [0, 1, 0],
        [0, 0, 1],
        [1, 0, 0]
    ])

    return pRt @ tRe @ eRd @ dRc @ cRb @ bRa @ aRw


def verifie_zone(x, y):
    return 45*x**2 + 30*x*y + 85*y**2 - 10.8*x - 8.4*y + 0.684 < 0


def dessiner_zone(lb, hg, hd, df_x, df_y, ti_x, ti_y):
    fig, ax = plt.subplots()

    # Contour de la pièce
    ax.plot([0, lb, lb, 0, 0], [hg, hg, 0, 0, hg], 'k-', lw=3)

    # Pente nominale de la tranche
    ax.plot([0, lb], [hg, hd], 'k--', lw=2)

    # Contour de la zone interdite
    x = np.linspace(-0.05, 0.2, 200)
    y = np.linspace(-0.05, 0.12, 200)
    X, Y = np.meshgrid(x, y)
    Z = 45*X**2 + 30*X*Y + 85*Y**2 - 10.8*X - 8.4*Y + 0.684
    ax.contour(X, Y, Z, [0], colors='red', lw=4)

    # Points de défauts
    ax.scatter(df_x, df_y, c='green', label='Défauts DF')

    # Points de mesure TI
    ax.scatter(ti_x, ti_y, c='blue', marker='x', s=60, label='Mesures TI')

    ax.set_xlabel('x (m)')
    ax.set_ylabel('y (m)')
    ax.grid(True)
    ax.set_aspect('equal')
    ax.legend()
    plt.title('Pièce avec zone interdite, défauts et mesures TI')
    plt.show()


def cinematique(theta_1, theta_2, theta_3, theta_4, theta_5, theta_6):
    """Calcule la cinématique directe du robot et retourne les positions de chaque articulation."""

    # Longueurs (mètres)
    l1, l2y, l3, l4y, l2x, l4x, l5, l6 = 0.15, 0.1, 0.5, 0.02, 0.05, 0.1, 0.3, 0.02

    # Vecteurs de position dans les repères w
    r_w_WwW = np.array([[0],   [0],   [0]])
    r_w_AwW = np.array([[0],   [l1],  [0]])
    r_a_BwA = np.array([[l2x], [l2y], [0]])
    r_b_CwB = np.array([[0],   [l3],  [0]])
    r_c_DwC = np.array([[l4x], [l4y], [0]])
    r_d_EwD = np.array([[l5],  [0],   [0]])
    r_e_TwE = np.array([[l6],  [0],   [0]])

    # matrices de rotation
    wRa, aRb, bRc, cRd, dRe, eRt = matrice_rotation(
        theta_1, theta_2, theta_3, theta_4, theta_5, theta_6
    )

    # Transformation des vecteurs dans le repère w
    r_w_BwA = wRa @ r_a_BwA
    r_w_CwB = wRa @ aRb @ r_b_CwB
    r_w_DwC = wRa @ aRb @ bRc @ r_c_DwC
    r_w_EwD = wRa @ aRb @ bRc @ cRd @ r_d_EwD
    r_w_TwE = wRa @ aRb @ bRc @ cRd @ dRe @ r_e_TwE

    # addition
    r_w_BwW = r_w_BwA + r_w_AwW
    r_w_CwW = r_w_CwB + r_w_BwA + r_w_AwW
    r_w_DwW = r_w_DwC + r_w_CwB + r_w_BwA + r_w_AwW
    r_w_EwW = r_w_EwD + r_w_DwC + r_w_CwB + r_w_BwA + r_w_AwW
    r_w_TwW = r_w_TwE + r_w_EwD + r_w_DwC + r_w_CwB + r_w_BwA + r_w_AwW

    matrice_final = np.hstack((r_w_WwW, r_w_AwW, r_w_BwW, r_w_CwW, r_w_DwW, r_w_EwW, r_w_TwW))

    return matrice_final, r_w_TwW


def dessiner_3D(matrice_final):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    Xs = matrice_final[0, :]
    Ys = matrice_final[1, :]
    Zs = matrice_final[2, :]

    ax.plot(Xs, Zs, Ys)

    ax.set_xlabel("X")
    ax.set_ylabel("z")
    ax.set_zlabel("y")
    ax.grid(True)
    ax.set_aspect('equal')
    plt.show()

def matrice_jacobienne(theta_1, theta_2, theta_3, theta_4, theta_5, theta_6):
    l1, l2y, l3, l4y, l2x, l4x, l5, l6 = 0.15, 0.1, 0.5, 0.02, 0.05, 0.1, 0.3, 0.02

    dxdth1 = -(l2x-l3*np.sin(theta_2)-l4y*np.sin(theta_2+theta_3)+(l4x+l5)*np.cos(theta_2+theta_3))*np.sin(theta_1)
    dzdth1 = -(l2x-l3*np.sin(theta_2)-l4y*np.sin(theta_2+theta_3)+(l4x+l5)*np.cos(theta_2+theta_3))*np.cos(theta_1)

    dxdth2 = -l3*np.cos(theta_2)-l4y*np.cos(theta_2+theta_3)-(l4x+l5)*np.sin(theta_2+theta_3)
    dydth2 = -l3*np.sin(theta_2)-l4y*np.sin(theta_2+theta_3)-(l4x+l5)*np.cos(theta_2+theta_3)

    dxdth3 = -l4y*np.cos(theta_2+theta_3)-(l4x+l5)*np.sin(theta_2+theta_3)
    dydth3 = -l4y*np.sin(theta_2+theta_3)-(l4x+l5)*np.cos(theta_2+theta_3)

    matrice_jacobienne = np.array([
        [dxdth1, dxdth2, dxdth3],
        [0, dydth2, dydth3],
        [dzdth1, 0, 0]
    ])

    matrice_jacobienne_A = np.array(
        [[0, dydth2, dydth3]])
    
    matrice_jacobienne_BC = np.array([
        [dxdth1, dxdth2, dxdth3],
        [0, dydth2, dydth3],
        [dzdth1, 0, 0]
    ])

    return matrice_jacobienne_A, matrice_jacobienne_BC

def cinematique_differentielle_A(theta_1, theta_2, theta_3, theta_4, theta_5, theta_6, vitesse):
    matrice_jacobienne_A, _ = matrice_jacobienne(theta_1, theta_2, theta_3, theta_4, theta_5, theta_6)
    matrice_vitesse = np.array([[vitesse]])

    # pas de l'inversion car c'est matrice de 1x1
    matrice_pseudo_inverse = matrice_jacobienne_A.T @ (matrice_jacobienne_A @ matrice_jacobienne_A.T)
    matrice_vitesse_angle = matrice_pseudo_inverse @ matrice_vitesse

    print(f"Situation 1: matrice de vitesse angulaire: \n{matrice_vitesse_angle}")


def cinematique_differentielle_B(theta_1, theta_2, theta_3, theta_4, theta_5, theta_6, vitesse):
    _, matrice_jacobienne_B = matrice_jacobienne(theta_1, theta_2, theta_3, theta_4, theta_5, theta_6)
    r = np.array([
        [0],
        [vitesse],
        [0]
    ])

    det_matrice_jacobienne_B = np.linalg.det(matrice_jacobienne_B)

    if det_matrice_jacobienne_B == 0:
        print("Erreur, determinant de matrice jacobienne égale à 0")
    else:
        matrice_vitesse_angle = np.linalg.inv(matrice_jacobienne_B) @ r
        print(f"Situation 2: matrice de vitesse angulaire: \n{matrice_vitesse_angle}")

def cinematique_differentielle_C(theta_1, theta_2, theta_3, theta_4, theta_5, theta_6, vitesse):
    _, matrice_jacobienne_C = matrice_jacobienne(theta_1, theta_2, theta_3, theta_4, theta_5, theta_6)

    matrice_JT_J = matrice_jacobienne_C.T @ matrice_jacobienne_C

    det_matrice_JT_J = np.linalg.det(matrice_JT_J)

    r = np.array([
        [0],
        [vitesse],
        [0]
    ])

    if det_matrice_JT_J == 0:
        print("Erreur, determinant de matrice jacobienne égale à 0")
    else:
        matrice_vitesse_angle = np.linalg.inv(matrice_JT_J) @ matrice_JT_J.T @ r
        round = np.round(matrice_vitesse_angle, 10)
        print(f"Situation 3: matrice de vitesse angulaire arrondie: \n{round}")
        print(f"Situation 3: matrice de vitesse angulaire non-arrondie: \n{matrice_vitesse_angle}")


def main():
    # Angles articulaires (rad) pour test cinematique
    # theta_1, theta_2, theta_3, theta_4, theta_5, theta_6 = 0.1, 0.1, 0.1, 0.1, 0.1, 0.1
    # theta_1, theta_2, theta_3, theta_4, theta_5, theta_6 = 0, 0, 0, 0, 0, 0

    # matrice_final, pt_t = cinematique(theta_1, theta_2, theta_3, theta_4, theta_5, theta_6)
    # print(f"Le point d'outil au frame w: {pt_t}")

    # dessiner_3D(matrice_final)
    # analyse_tranche()

    # Angles articulaires (rad) pour test cinematique diff
    theta_1, theta_2, theta_3, theta_4, theta_5, theta_6 = -0.4, -1.2, 0, 0, -0.3708, 0
    # theta_1, theta_2, theta_3, theta_4, theta_5, theta_6 = 0, 0, 1.521, 0, 0, 0
    vitesse = 1
    cinematique_differentielle_A(theta_1, theta_2, theta_3, theta_4, theta_5, theta_6, vitesse)
    cinematique_differentielle_B(theta_1, theta_2, theta_3, theta_4, theta_5, theta_6, vitesse)
    cinematique_differentielle_C(theta_1, theta_2, theta_3, theta_4, theta_5, theta_6, vitesse)


main()
