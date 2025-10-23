# GemmaGonzalez142231_Trabajo2.R
# Trabajo final Bioinformática - Curso 25/26
# Análisis de parámetros biomédicos por tratamiento

# 1. Cargar librerías (si necesarias) y datos del archivo "datos_biomed.csv". (0.5 pts)

#cargamos la lireria en la que tenemos ggplot para hacer las diferentes gráficas 
library(ggplot2)

#y ahora leemos el archivo y lo guardamos en una variable para que sea más facil trabajar con ella
datos <- read.csv("datos_biomed.csv")

# 2. Exploración inicial con las funciones head(), summary(), dim() y str(). ¿Cuántas variables hay? ¿Cuántos tratamientos? (0.5 pts)

#head sirve para ver las primeras 6 filas
head(datos)

#summary para ver el resumen estadístico de los datos
summary(datos)

#dim te muestra las filas y columnas del archivo
dim(datos)

#str te muestra los tipos de datos que se encuentran en el archivo
str(datos)

#para ver el numero de variables iriamos directamente al numero de columnas 
n_variables <- ncol(datos)

#para ver el numero de tratamientos que hay vamos directamente a los datos de la columna de tratamientos, y para que no haya repeticiones ponemos unique
n_tratamiento <- length(unique(datos$Tratamiento))

#y llamamos a las funciones para ver los resultados
n_variables
n_tratamiento

# 3. Una gráfica que incluya todos los boxplots por tratamiento. (1 pt)
#creamos el grafico, indicando que el tratamiento vaya en el eje x, y los niveles de glucosa en el eje y, y fill nos indica el color, en este caso de como se ve el tratamiento en la grafica
#con geom_boxplot indico que sea un grafico de boxplot
#con theme minimal indico que solo quiero que sea visual
ggplot(datos, aes(x = Tratamiento, y = Glucosa, fill = Tratamiento)) +
	geom_boxplot() + 
	theme_minimal() +
	labs(title = "Concentración de Glucosa por tratamiento")

# 4. Realiza un violin plot (investiga qué es). (1 pt)
#este gráfico se trata de una mezcla de boxplot y una curva de densidad, es decir, boxplot porque muestra mediana y cuartiles, y de forma de violin porque con la curva de densidad indica que las zonas mas estrechas tienen menor concentración de datos, y las más anchas más concetración.
#para hacerlo de tipo violin pondremos geom_violin, pero dentro de el indicaremos trim False para que no se recorten las colas
#width es la anchura
ggplot(datos, aes(x = Tratamiento, y = Glucosa, fill = Tratamiento)) +
	geom_violin(trim = FALSE) +
	theme_minimal() +
	labs(theme = "Concentración de Glucosa por Tratamiento")

# 5. Realiza un gráfico de dispersión "Glucosa vs Presión". Emplea legend() para incluir una leyenda en la parte inferior derecha. (1 pt)
#hacemos un gráfico de dispersion simple con los datos de glucosa y presion
#ponemos la glucosa en el eje x y la presión en el eje y, col seria para poner un color distinto a tratamiento, con pch indicas como quieres que se vean los puntos, y he elegido 19 para que sean redonditos
plot(datos$Glucosa, datos$Presion,
	col = as.factor(datos$Tratamiento),
	pch = 19,
	main = "Glucosa vs Presión por Tratamiento",
	xlab = "Glucosa",
	ylab = "Presión")

#Y ahora agregamos una leyenda, y lo añadimos abajo en la gráfica
legend("bottomright",
	legend = unique(datos$Tratamiento),
	col = 1:length(unique(datos$Tratamiento)),
	pch = 19,
	title = "Tratamiento")

# 6. Realiza un facet Grid (investiga qué es): Colesterol vs Presión por tratamiento. (1 pt)
#lo que hace un face grid es dividir el gráfico por paneles, por columnas y filas. En este caso se ven en los paneles los diferentes tratamientos
ggplot(datos, aes(x = Presion, y = Colesterol)) +
	geom_point(aes(color = Tratamiento)) +
	facet_grid(~Tratamiento) +
	theme_minimal() +
	labs(title = "Colesterol vs Presión por Tratamiento")

# 7. Realiza un histogramas para cada variable. (0.5 pts)
#dodge hace que las barras no se superpongan entre sí, y bins divide los datos en intervalos, yo he puesto 20 porque asi se vene mas barras y es mas detallado
ggplot(datos, aes(x = Glucosa, fill = Tratamiento)) +
	geom_histogram(position = "dodge", bins = 20, color = "black") +
	theme_minimal() +
	labs(title = "Glucosa por Tratamiento")

ggplot(datos, aes(x = Presion, fill = Tratamiento)) +
	geom_histogram(position = "dodge", bins = 20, color = "black") +
	theme_minimal() +
	labs(title = "Presion por Tratamiento")

ggplot(datos, aes(x = Colesterol, fill = Tratamiento)) +
	geom_histogram(position = "dodge", bins = 20, color = "black") +
	theme_minimal() +
	labs(title = "Colesterol por Tratamiento")


# 8. Crea un factor a partir del tratamiento. Investifa factor(). (1 pt)
#con factor se convierte una variable, en este caso numerica, a factor, para que sea más fácil de manejar para funciones estadísticas, ya que lo identifican como grupos
datos$Tratamiento <- factor(datos$Tratamiento)

# 9. Obtén la media y desviación estándar de los niveles de glucosa por tratamiento. Emplea aggregate() o apply(). (0.5 pts)
#a parte de aggregate, con fun indicamos que funciones queremos ver, como media y desviación
datos_glucosa <- aggregate(Glucosa ~ Tratamiento, data = datos,
	FUN = function(x) c(media = mean(x), sd = sd(x)))

#y lo llamamos para ver los resultados
datos_glucosa 

# 10. Extrae los datos para cada tratamiento y almacenalos en una variable. Ejemplo todos los datos de Placebo en una variable llamada placebo. (1 pt)
#para extraerlos podemos llamarlos directamente y definirles una variable, extraigo de los datos la columna tratamiento, y de ahí extraigo los datos asociados al nombre Placebo. Y asi con todos los fármacos.
placebo <- datos[datos$Tratamiento == "Placebo", ]
tratamientoA <- datos[datos$Tratamiento == "FarmacoA", ]
tratamientoB <- datos[datos$Tratamiento == "FarmacoB", ]

#y los llamamos para ver sus resultados
placebo
tratamientoA
tratamientoB

# 11. Evalúa si los datos siguen una distribución normal y realiza una comparativa de medias acorde. (1 pt)
#para medir la distribucion podemos utilizar el shapiro test, que puede aplicar la distribución por cada tratamiento
#tambien podremos ver el pvalue, por lo que si es menor de 0.05, en este caso rechazamos la hipotesis nula
by(datos$Glucosa, datos$Tratamiento, shapiro.test)
by(datos$Presion, datos$Tratamiento, shapiro.test)
by(datos$Colesterol, datos$Tratamiento, shapiro.test)

#se puede ver que todos los pvalues superan el 0.05, por lo que todos muestran una distribución normal
#en este caso se cumple la hipotesis nula de que siguen una distribución normal, y además en todos los casos
#si cumple la hipotesis nula podemos utilizar t.student con t , para realizar la comparativa
t.test(placebo$Glucosa, tratamientoA$Glucosa)

t.test(placebo$Glucosa, tratamientoB$Glucosa)

t.test(tratamientoA$Glucosa, tratamientoB$Glucosa)

t.test(placebo$Presion, tratamientoA$Presion)

t.test(placebo$Presion, tratamientoB$Presion)

t.test(tratamientoA$Presion, tratamientoB$Presion)

t.test(placebo$Colesterol, tratamientoA$Colesterol)

t.test(placebo$Colesterol, tratamientoB$Colesterol)

t.test(tratamientoA$Colesterol, tratamientoB$Colesterol)


#como conclusiones al respecto, en el caso de la glucosa, su p-value de comparacion de medias muestra un valor mayor a 0,05, por lo que se aprueba la hipotesis nula, que todas sus medias son similares entre los grupos
#en el caso de la presion y el colesterol, si hay diferencias significativas ya que ambos tienen una pvalue menor que 0.05, por lo que se rechaza la hipotesis nula, y se confirma que sus medias son diferentes en cada grupo
#(menos en los casos de comparación entre el tratamiento A respecto al B teniendo en cuenta la Presión, esta tiene medias similares, y comparación del placebo con el tratamiento B en cuanto al colesterol.

# 12. Realiza un ANOVA sobre la glucosa para cada tratamiento. (1 pt)
#es igual que el ejercicio anterior, pero aplicando anova, comparando medias
anova_glucosa <- aov(Glucosa ~ Tratamiento, data = datos)
summary(anova_glucosa)

#el Pr(>F)indica el pvalue, y como es mayor a 0.05, podemos decir que las medias de los grupos son similares


