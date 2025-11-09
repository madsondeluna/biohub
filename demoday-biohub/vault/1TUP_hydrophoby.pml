# Script PyMOL gerado pelo BioHub
# Propriedade: hydrophobicity

# Carrega a estrutura
load /Users/madsonluna/Downloads/biohub-alpha/1TUP_hydrophoby.pdb, 1TUP_hydrophoby

# Remove todas as representações padrão
hide everything, 1TUP_hydrophoby

# === HIDROFOBICIDADE ===
# Azul = Hidrofílico (-4.5), Vermelho = Hidrofóbico (4.5)

# Representação Cartoon (fita)
show cartoon, 1TUP_hydrophoby
cartoon automatic, 1TUP_hydrophoby
set cartoon_fancy_helices, 1
spectrum b, blue_white_red, 1TUP_hydrophoby, minimum=-4.5, maximum=4.5

# Representação Sticks (bastões)
show sticks, 1TUP_hydrophoby
set stick_radius, 0.2, 1TUP_hydrophoby
set stick_color, gray, 1TUP_hydrophoby
util.cbag 1TUP_hydrophoby

# Representação Surface (superfície)
show surface, 1TUP_hydrophoby
set surface_quality, 1
set transparency, 0.5, 1TUP_hydrophoby
# Aplica gradiente de hidrofobicidade na superfície
set surface_color, white, 1TUP_hydrophoby
spectrum b, blue_white_red, 1TUP_hydrophoby, minimum=-4.5, maximum=4.5

# === Configurações Gerais de Qualidade ===
bg_color white
set ray_shadow, 0
set antialias, 2
set orthoscopic, 0
set valence, 0

# Centra e ajusta visualização
center 1TUP_hydrophoby
zoom 1TUP_hydrophoby
orient 1TUP_hydrophoby

# Comandos úteis:
# hide surface - ocultar superfície
# hide sticks - ocultar bastões
# hide cartoon - ocultar cartoon
# set transparency, 0.7 - aumentar transparência

# Para salvar a sessão, use:
# save /Users/madsonluna/Downloads/biohub-alpha/1TUP_hydrophoby.pse