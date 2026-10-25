"""
Actividad formativa 2: Simulacion y analisis de senales con Fourier.

Generacion de senales y representacion en el tiempo.
Calculo y visualizacion de magnitud y fase mediante FFT.
Verificacion de linealidad, desplazamiento temporal y escalamiento.
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def crear_senales():
    """Genera las tres senales sobre un mismo eje de tiempo."""

    # Tomamos 1000 muestras por segundo durante 2 segundos.
    frecuencia_muestreo = 1000  # Hz
    duracion = 2.0             # segundos
    numero_muestras = int(frecuencia_muestreo * duracion)

    # Eje temporal desde -1 s hasta 0.999 s, con pasos de 0.001 s.
    tiempo = (
                     np.arange(numero_muestras) - numero_muestras // 2
             ) / frecuencia_muestreo

    # Pulso de amplitud 1, centrado en cero y de ancho 0.4 s.
    ancho_pulso = 0.4
    pulso = np.where(np.abs(tiempo) < ancho_pulso / 2, 1.0, 0.0)

    # Escalon: vale 0 antes de t = 0 y 1 desde t = 0.
    escalon = np.where(tiempo >= 0, 1.0, 0.0)

    # Senoidal de amplitud 1 y frecuencia de 5 Hz.
    frecuencia_senoidal = 5.0
    senoidal = np.sin(2 * np.pi * frecuencia_senoidal * tiempo)

    senales = {
        "Pulso rectangular: ancho de 0.4 s": pulso,
        "Escalon unitario": escalon,
        "Senal senoidal: frecuencia de 5 Hz": senoidal,
    }

    return tiempo, senales


def graficar_senales(tiempo, senales):
    """Muestra las senales y guarda una imagen con las graficas."""

    figura, ejes = plt.subplots(
        nrows=3,
        ncols=1,
        figsize=(10, 8),
        sharex=True,
        constrained_layout=True,
    )

    colores = ["tab:blue", "tab:orange", "tab:green"]

    for eje, (nombre, valores), color in zip(
            ejes, senales.items(), colores
    ):
        eje.plot(tiempo, valores, color=color, linewidth=1.8)
        eje.set_title(nombre)
        eje.set_ylabel("Amplitud")
        eje.grid(True, alpha=0.3)
        eje.axhline(0, color="gray", linewidth=0.7)
        eje.axvline(0, color="gray", linewidth=0.7)
        eje.set_xlim(tiempo[0], tiempo[-1])

    ejes[-1].set_xlabel("Tiempo (s)")
    figura.suptitle("Senales en el dominio del tiempo", fontsize=16)

    # La carpeta se crea junto a main.py, si todavia no existe.
    carpeta_resultados = Path(__file__).resolve().parent / "resultados"
    carpeta_resultados.mkdir(exist_ok=True)

    archivo_imagen = carpeta_resultados / "senales_tiempo.png"
    figura.savefig(archivo_imagen, dpi=200)

    print(f"Grafica guardada en: {archivo_imagen}")
    plt.show()

def calcular_espectro(tiempo, valores):
    """Aproxima la transformada de Fourier de la senal observada."""

    intervalo_muestreo = tiempo[1] - tiempo[0]
    numero_muestras = len(valores)

    # Frecuencias correspondientes a los resultados de la FFT.
    frecuencias = np.fft.fftshift(
        np.fft.fftfreq(numero_muestras, d=intervalo_muestreo)
    )

    # Nuestro eje temporal tiene t = 0 en el centro del arreglo.
    # ifftshift coloca esa muestra al inicio antes de calcular la FFT.
    valores_ordenados = np.fft.ifftshift(valores)

    # Multiplicar por el intervalo de muestreo aproxima la integral
    # de Fourier de la senal dentro del intervalo observado.
    espectro = intervalo_muestreo * np.fft.fftshift(
        np.fft.fft(valores_ordenados)
    )

    return frecuencias, espectro


def graficar_espectros(tiempo, senales):
    """Grafica la magnitud y la fase del espectro de cada senal."""

    figura, ejes = plt.subplots(
        nrows=3,
        ncols=2,
        figsize=(12, 9),
        constrained_layout=True,
    )

    colores = ["tab:blue", "tab:orange", "tab:green"]

    for fila, ((nombre, valores), color) in enumerate(
            zip(senales.items(), colores)
    ):
        frecuencias, espectro = calcular_espectro(tiempo, valores)

        magnitud = np.abs(espectro)
        fase = np.angle(espectro)

        # La fase no esta definida cuando la magnitud es cero.
        # Ocultamos tambien los valores cercanos al ruido numerico.
        umbral = magnitud.max() * 1e-6
        fase = np.where(magnitud > umbral, fase, np.nan)

        # Calculamos todo el espectro, pero mostramos de -20 a 20 Hz
        # para apreciar mejor sus componentes principales.
        zona_visible = np.abs(frecuencias) <= 20

        eje_magnitud = ejes[fila, 0]
        eje_fase = ejes[fila, 1]

        eje_magnitud.plot(
            frecuencias[zona_visible],
            magnitud[zona_visible],
            color=color,
            marker=".",
            markersize=4,
            linewidth=1.2,
        )
        eje_magnitud.set_title(f"{nombre}\nMagnitud", fontsize=10)
        eje_magnitud.set_ylabel("|X(f)| (amplitud x s)")
        eje_magnitud.set_ylim(bottom=0)

        # Usamos puntos para no unir saltos de fase ni zonas indefinidas.
        eje_fase.plot(
            frecuencias[zona_visible],
            fase[zona_visible],
            linestyle="none",
            marker=".",
            markersize=5,
            color=color,
        )
        eje_fase.set_title(f"{nombre}\nFase", fontsize=10)
        eje_fase.set_ylabel("Fase (rad)")
        eje_fase.set_yticks(
            [-np.pi, -np.pi / 2, 0, np.pi / 2, np.pi],
            ["-pi", "-pi/2", "0", "pi/2", "pi"],
        )
        eje_fase.set_ylim(-np.pi - 0.3, np.pi + 0.3)

        for eje in (eje_magnitud, eje_fase):
            eje.set_xlabel("Frecuencia (Hz)")
            eje.set_xlim(-20, 20)
            eje.grid(True, alpha=0.3)
            eje.axvline(0, color="gray", linewidth=0.7)

    figura.suptitle(
        "Espectros de las senales observadas: magnitud y fase",
        fontsize=15,
    )

    carpeta_resultados = Path(__file__).resolve().parent / "resultados"
    carpeta_resultados.mkdir(exist_ok=True)

    archivo_imagen = carpeta_resultados / "espectros_fourier.png"
    figura.savefig(archivo_imagen, dpi=200)

    print(f"Espectros guardados en: {archivo_imagen}")
    plt.show()

def verificar_linealidad(tiempo, senales):
    """Compara la FFT de una combinacion con la combinacion de FFT."""

    pulso = senales["Pulso rectangular: ancho de 0.4 s"]
    senoidal = senales["Senal senoidal: frecuencia de 5 Hz"]

    a = 2.0
    b = 0.5

    # Primera ruta: combinar las senales y transformar el resultado.
    combinacion = a * pulso + b * senoidal
    frecuencias, espectro_combinacion = calcular_espectro(
        tiempo, combinacion
    )

    # Segunda ruta: transformar por separado y combinar los espectros.
    _, espectro_pulso = calcular_espectro(tiempo, pulso)
    _, espectro_senoidal = calcular_espectro(tiempo, senoidal)

    espectro_referencia = (
            a * espectro_pulso + b * espectro_senoidal
    )

    # Comparamos los espectros complejos, incluyendo magnitud y fase.
    diferencia = espectro_combinacion - espectro_referencia
    error_relativo = (
            np.linalg.norm(diferencia)
            / np.linalg.norm(espectro_referencia)
    )

    cumple = np.allclose(
        espectro_combinacion,
        espectro_referencia,
        rtol=1e-10,
        atol=1e-12,
    )

    print("\nCOMPROBACION DE LINEALIDAD")
    print(f"Error relativo: {error_relativo:.3e}")
    print(f"Coinciden dentro de la tolerancia: {cumple}")

    figura, ejes = plt.subplots(
        nrows=2,
        ncols=1,
        figsize=(10, 8),
        constrained_layout=True,
    )

    # Representacion temporal de la combinacion.
    ejes[0].plot(tiempo, a * pulso, "--", label="2 x pulso")
    ejes[0].plot(tiempo, b * senoidal, ":", label="0.5 x senoidal")
    ejes[0].plot(tiempo, combinacion, label="Combinacion", linewidth=1.8)
    ejes[0].set_title("Combinacion de las senales en el tiempo")
    ejes[0].set_xlabel("Tiempo (s)")
    ejes[0].set_ylabel("Amplitud")

    # Superponemos las magnitudes obtenidas por las dos rutas.
    zona_visible = np.abs(frecuencias) <= 20

    ejes[1].plot(
        frecuencias[zona_visible],
        np.abs(espectro_combinacion[zona_visible]),
        label="Transformada de la combinacion",
        linewidth=1.8,
    )
    ejes[1].plot(
        frecuencias[zona_visible],
        np.abs(espectro_referencia[zona_visible]),
        linestyle="none",
        marker="o",
        fillstyle="none",
        markersize=4,
        label="Combinacion de las transformadas",
    )
    ejes[1].set_title("Comparacion de los espectros")
    ejes[1].set_xlabel("Frecuencia (Hz)")
    ejes[1].set_ylabel("|X(f)| (amplitud x s)")
    ejes[1].set_xlim(-20, 20)

    for eje in ejes:
        eje.grid(True, alpha=0.3)
        eje.legend()

    figura.suptitle("Propiedad de linealidad", fontsize=16)

    carpeta_resultados = Path(__file__).resolve().parent / "resultados"
    carpeta_resultados.mkdir(exist_ok=True)

    archivo_imagen = carpeta_resultados / "linealidad.png"
    figura.savefig(archivo_imagen, dpi=200)

    print(f"Grafica guardada en: {archivo_imagen}")
    plt.show()


def verificar_desplazamiento(tiempo, senales):
    """Comprueba el efecto de retrasar una senal en el tiempo."""

    pulso = senales["Pulso rectangular: ancho de 0.4 s"]

    intervalo_muestreo = tiempo[1] - tiempo[0]
    retraso_deseado = 0.25  # segundos

    # Convertimos el retraso a un numero entero de muestras.
    muestras_retraso = int(round(
        retraso_deseado / intervalo_muestreo
    ))
    retraso_real = muestras_retraso * intervalo_muestreo

    # Verificamos que el desplazamiento no corte parte del pulso.
    if np.any(pulso[-muestras_retraso:] != 0):
        raise ValueError("El pulso desplazado saldria del intervalo.")

    # Desplazamos hacia la derecha y rellenamos el inicio con ceros.
    pulso_retrasado = np.zeros_like(pulso)
    pulso_retrasado[muestras_retraso:] = pulso[:-muestras_retraso]

    frecuencias, espectro_original = calcular_espectro(tiempo, pulso)
    _, espectro_retrasado = calcular_espectro(tiempo, pulso_retrasado)

    # Propiedad: x(t - t0) corresponde a X(f) * exp(-j*2*pi*f*t0).
    factor_fase = np.exp(-2j * np.pi * frecuencias * retraso_real)
    espectro_teorico = espectro_original * factor_fase

    error_relativo = (
            np.linalg.norm(espectro_retrasado - espectro_teorico)
            / np.linalg.norm(espectro_teorico)
    )

    cumple = np.allclose(
        espectro_retrasado,
        espectro_teorico,
        rtol=1e-10,
        atol=1e-12,
    )

    magnitud_conservada = np.allclose(
        np.abs(espectro_retrasado),
        np.abs(espectro_original),
        rtol=1e-10,
        atol=1e-12,
    )

    print("\nCOMPROBACION DE DESPLAZAMIENTO TEMPORAL")
    print(f"Retraso aplicado: {retraso_real:.3f} segundos")
    print(f"Muestras de retraso: {muestras_retraso}")
    print(f"Error relativo: {error_relativo:.3e}")
    print(f"Coincide con la teoria: {cumple}")
    print(f"Magnitud conservada: {magnitud_conservada}")

    figura, ejes = plt.subplots(
        nrows=3,
        ncols=1,
        figsize=(10, 9),
        constrained_layout=True,
    )

    # Primera grafica: posicion original y posicion retrasada.
    ejes[0].plot(tiempo, pulso, label="Pulso original")
    ejes[0].plot(
        tiempo, pulso_retrasado, "--", label="Pulso retrasado 0.25 s"
    )
    ejes[0].set_title("Desplazamiento del pulso en el tiempo")
    ejes[0].set_xlabel("Tiempo (s)")
    ejes[0].set_ylabel("Amplitud")

    # Segunda grafica: las magnitudes deben coincidir.
    zona_visible = np.abs(frecuencias) <= 20

    ejes[1].plot(
        frecuencias[zona_visible],
        np.abs(espectro_original[zona_visible]),
        label="Magnitud original",
    )
    ejes[1].plot(
        frecuencias[zona_visible],
        np.abs(espectro_retrasado[zona_visible]),
        linestyle="none",
        marker="o",
        fillstyle="none",
        markersize=4,
        label="Magnitud del pulso retrasado",
    )
    ejes[1].set_title("La magnitud del espectro se conserva")
    ejes[1].set_xlabel("Frecuencia (Hz)")
    ejes[1].set_ylabel("|X(f)| (amplitud x s)")
    ejes[1].set_xlim(-20, 20)

    # Tercera grafica: cambio de fase introducido por el retraso.
    # Evitamos dividir entre valores espectrales casi nulos.
    umbral = np.max(np.abs(espectro_original)) * 1e-6
    puntos_validos = zona_visible & (
            np.abs(espectro_original) > umbral
    )
    frecuencias_fase = frecuencias[puntos_validos]

    cociente = (
            espectro_retrasado[puntos_validos]
            / espectro_original[puntos_validos]
    )

    # Eliminamos saltos de 2*pi para mostrar la pendiente de la fase.
    cambio_fase = np.unwrap(np.angle(cociente))
    indice_cero = np.argmin(np.abs(frecuencias_fase))
    cambio_fase = cambio_fase - cambio_fase[indice_cero]

    fase_teorica = -2 * np.pi * frecuencias_fase * retraso_real

    ejes[2].plot(
        frecuencias_fase,
        fase_teorica,
        label="Cambio de fase teorico: -2*pi*f*t0",
    )
    ejes[2].plot(
        frecuencias_fase,
        cambio_fase,
        linestyle="none",
        marker="o",
        fillstyle="none",
        markersize=4,
        label="Cambio de fase calculado",
    )
    ejes[2].set_title("Cambio de fase sin saltos de 2*pi")
    ejes[2].set_xlabel("Frecuencia (Hz)")
    ejes[2].set_ylabel("Cambio de fase (rad)")
    ejes[2].set_xlim(-20, 20)

    for eje in ejes:
        eje.grid(True, alpha=0.3)
        eje.legend()

    figura.suptitle("Propiedad de desplazamiento temporal", fontsize=15)

    carpeta_resultados = Path(__file__).resolve().parent / "resultados"
    carpeta_resultados.mkdir(exist_ok=True)

    archivo_imagen = carpeta_resultados / "desplazamiento_temporal.png"
    figura.savefig(archivo_imagen, dpi=200)

    print(f"Grafica guardada en: {archivo_imagen}")
    plt.show()


def verificar_escalamiento(tiempo):
    """Comprueba el escalamiento usando una senal gaussiana."""

    parametro_ancho = 0.2  # segundos; controla el ancho de la campana
    factor = 2.0

    # Senal original: x(t) = exp(-pi * (t / ancho)^2).
    original = np.exp(-np.pi * (tiempo / parametro_ancho) ** 2)

    # Senal comprimida: y(t) = x(2t).
    comprimida = np.exp(
        -np.pi * (factor * tiempo / parametro_ancho) ** 2
    )

    frecuencias, espectro_original = calcular_espectro(tiempo, original)
    _, espectro_comprimido = calcular_espectro(tiempo, comprimida)

    def transformada_gaussiana(f):
        """Formula analitica de la transformada de la senal original."""
        return parametro_ancho * np.exp(
            -np.pi * (parametro_ancho * f) ** 2
        )

    # Referencia analitica para comprobar la FFT original.
    referencia_original = transformada_gaussiana(frecuencias)

    # Aplicamos la propiedad: Y(f) = X(f / a) / |a|.
    referencia_comprimida = (
            transformada_gaussiana(frecuencias / factor) / abs(factor)
    )

    error_original = (
            np.linalg.norm(espectro_original - referencia_original)
            / np.linalg.norm(referencia_original)
    )
    error_escalamiento = (
            np.linalg.norm(espectro_comprimido - referencia_comprimida)
            / np.linalg.norm(referencia_comprimida)
    )

    cumple_original = np.allclose(
        espectro_original,
        referencia_original,
        rtol=1e-10,
        atol=1e-12,
    )
    cumple_escalamiento = np.allclose(
        espectro_comprimido,
        referencia_comprimida,
        rtol=1e-10,
        atol=1e-12,
    )

    print("\nCOMPROBACION DE ESCALAMIENTO")
    print(f"Factor aplicado: {factor:.1f}")
    print(f"Error relativo de la original: {error_original:.3e}")
    print(f"Error relativo de escalamiento: {error_escalamiento:.3e}")
    print(f"Original coincide con la teoria: {cumple_original}")
    print(f"Escalamiento coincide con la teoria: {cumple_escalamiento}")

    figura, ejes = plt.subplots(
        nrows=2,
        ncols=1,
        figsize=(10, 8),
        constrained_layout=True,
    )

    # Comparacion temporal: la segunda campana es mas estrecha.
    ejes[0].plot(tiempo, original, label="Original: x(t)")
    ejes[0].plot(tiempo, comprimida, "--", label="Comprimida: x(2t)")
    ejes[0].set_title("Compresion de la senal en el tiempo")
    ejes[0].set_xlabel("Tiempo (s)")
    ejes[0].set_ylabel("Amplitud")
    ejes[0].set_xlim(-0.5, 0.5)

    # Comparacion espectral: la senal comprimida tiene mayor ancho
    # en frecuencia y menor magnitud maxima.
    zona_visible = np.abs(frecuencias) <= 15

    ejes[1].plot(
        frecuencias[zona_visible],
        np.abs(espectro_original[zona_visible]),
        label="Espectro original: |X(f)|",
    )
    ejes[1].plot(
        frecuencias[zona_visible],
        np.abs(espectro_comprimido[zona_visible]),
        "--",
        label="Espectro calculado de x(2t)",
    )
    ejes[1].plot(
        frecuencias[zona_visible],
        np.abs(referencia_comprimida[zona_visible]),
        linestyle="none",
        marker="o",
        fillstyle="none",
        markersize=4,
        label="Prediccion teorica: |X(f/2)| / 2",
    )
    ejes[1].set_title("Expansion del espectro en frecuencia")
    ejes[1].set_xlabel("Frecuencia (Hz)")
    ejes[1].set_ylabel("Magnitud (amplitud x s)")
    ejes[1].set_xlim(-15, 15)

    for eje in ejes:
        eje.grid(True, alpha=0.3)
        eje.legend()

    figura.suptitle("Propiedad de escalamiento", fontsize=16)

    carpeta_resultados = Path(__file__).resolve().parent / "resultados"
    carpeta_resultados.mkdir(exist_ok=True)

    archivo_imagen = carpeta_resultados / "escalamiento.png"
    figura.savefig(archivo_imagen, dpi=200)

    print(f"Grafica guardada en: {archivo_imagen}")
    plt.show()


def main():
    tiempo, senales = crear_senales()
    graficar_senales(tiempo, senales)
    graficar_espectros(tiempo, senales)
    verificar_linealidad(tiempo, senales)
    verificar_desplazamiento(tiempo, senales)
    verificar_escalamiento(tiempo)


if __name__ == "__main__":
    main()