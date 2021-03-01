library("sessioninfo")
library("here")
library("ggplot2")

## Hello world
print("Soy Leo")

## Crear directorio para las figuras
dir.create(here::here("figuras"), showWarnings = FALSE)

## Hacer una imagen de ejemplo
pdf(here::here("figuras", "mtcars_gear_vs_mpg.pdf"),
    useDingbats = FALSE
)
ggplot(mtcars, aes(group = gear, y = mpg)) +
  geom_boxplot()
dev.off()

## Para reproducir mi código
options(width = 120)
sessioninfo::session_info()
