# Actividad-2-SCII
Códigos correspondientes al desarrollo del informe 2 de la materia Sistemas de control II, control por realimentación de estados
En este archivo se detalla la funcionalidad de cada código como orientación:

MOTOR1: Primer intento de control del motor, no termina de generar una respuesta sin error en estado estacionario
MOTOR1_masintegrador: Se agrega un integrador al sistema, se produce una buena respuesta, todavía no hay observador
MOTOR2_observador: Código que se usó para pruebas, ignorar, pasar a "MOTOR2_compar_con_obs"
MOTOR2_compar_con_obs: Se agrega el observador cuando no puede medirse la corriente y se compara la respuesta con el sistema realimentado con las variables reales

AVION1: Control del sistema de cabeceo del avión con polos asignados, sin observador
AVION2_observador: Se agrega un observador midiendo solamente la altura mas la comparación con el sistema realimentado con las variables reales

PENDULO_1: Primer control del pendulo sin gráficos muy detallados, mirar "PENDULO_2_OBSERVADOR"
PENDULO_2_OBSERVADOR: Simulación del pendulo en el equilibrio inestable del sistema lineal con y sin observador para comparación. Valores de Q definidos por iteracion en "iteracion_controladores_pendulo2"
iteracion_controladores_pendulo2: Iteración de controladores K calculados de manera aleatoria (previo conociendo aproximadamente las magnitudes de Q) para encontrar por "fuerza bruta" un controlador óptimo que no produzca mucha oscilación de tita
PENDULO_2_OBSERVADOR_nl: Aplicación del controlador y observador utilizados en "PENDULO_2_OBSERVADOR" en el sistema no lineal
PENDULO_3: Primer control del pendulo sin gráficos muy detallados y sin integrador, código inicial del cual partieron las soluciones de la grúa
PENDULO_3_integrador: Evolución del código "PENDULO_3" agregando integrador, aún se puede mejorar con el cambio de controladores
PENDULO_3_integrador_dos_controladores: Prueba de switcheo de controladores para cambios de masa, mejora respecto a "PENDULO_3_integrador", aunque se sigue buscando el controlador más adecuado para cuando se agregue el observador:
PENDULO_4_integrador_observador: Versión final del control de la grúa linealizada. Valores de Q definidos por iteracion en "iteracion_controladores_pendulo4"
iteracion_controladores_pendulo4: Iteración de controladores K calculados de manera aleatoria (previo conociendo aproximadamente las magnitudes de Q) para encontrar por "fuerza bruta" un controlador óptimo que no produzca mucha oscilación de tita aún cuando se observa al sistema (osc. max <0.05rad)
PENDULO_4_integrador_observador_nl: Aplicación de los controladores y observador utilizados en "PENDULO_4_integrador_observador" sobre el sistema no lineal de la grúa
