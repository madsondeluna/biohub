# Script PyMOL gerado pelo BioHub
# Propriedade: hydrophobicity

# Carrega a estrutura
load /Users/madsonluna/Documents/biohub-1/4HHB_hydro.pdb, 4HHB_hydro

# Remove todas as representações padrão
hide everything, 4HHB_hydro

# === HIDROFOBICIDADE ===
# Azul = Hidrofílico (-4.5), Vermelho = Hidrofóbico (4.5)

# Representação Cartoon (fita)
show cartoon, 4HHB_hydro
cartoon automatic, 4HHB_hydro
set cartoon_fancy_helices, 1
spectrum b, blue_white_red, 4HHB_hydro, minimum=-4.5, maximum=4.5

# Representação Sticks (bastões)
show sticks, 4HHB_hydro
set stick_radius, 0.2, 4HHB_hydro
set stick_color, gray, 4HHB_hydro
util.cbag 4HHB_hydro

# Representação Surface (superfície)
show surface, 4HHB_hydro
set surface_quality, 1
set transparency, 0.5, 4HHB_hydro
# Aplica gradiente de hidrofobicidade na superfície
set surface_color, white, 4HHB_hydro
spectrum b, blue_white_red, 4HHB_hydro, minimum=-4.5, maximum=4.5

# === Configurações Gerais de Qualidade ===
bg_color white
set ray_shadow, 0
set antialias, 2
set orthoscopic, 0
set valence, 0

# Centra e ajusta visualização
center 4HHB_hydro
zoom 4HHB_hydro
orient 4HHB_hydro

# Comandos úteis:
# hide surface - ocultar superfície
# hide sticks - ocultar bastões
# hide cartoon - ocultar cartoon
# set transparency, 0.7 - aumentar transparência

# Para salvar a sessão, use:
# save /Users/madsonluna/Documents/biohub-1/4HHB_hydro.pse