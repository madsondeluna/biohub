# Script PyMOL gerado pelo BioHub
# Propriedade: sasa

# Carrega a estrutura
load /Users/madsonluna/Downloads/biohub-alpha/1TUP_sasa.pdb, 1TUP_sasa

# Remove todas as representações padrão
hide everything, 1TUP_sasa

# === SASA (Acessibilidade ao Solvente) ===
# Vermelho = Enterrado (0.0 Ų), Azul = Exposto (10.2 Ų)

# Representação Cartoon (fita)
show cartoon, 1TUP_sasa
cartoon automatic, 1TUP_sasa
set cartoon_fancy_helices, 1
spectrum b, red_white_blue, 1TUP_sasa, minimum=0.0, maximum=10.240186816308773

# Representação Sticks (bastões)
show sticks, 1TUP_sasa
set stick_radius, 0.2, 1TUP_sasa
set stick_color, gray, 1TUP_sasa
util.cbag 1TUP_sasa

# Representação Surface (superfície)
show surface, 1TUP_sasa
set surface_quality, 1
set transparency, 0.5, 1TUP_sasa
# Aplica gradiente de SASA na superfície
set surface_color, white, 1TUP_sasa
spectrum b, red_white_blue, 1TUP_sasa, minimum=0.0, maximum=10.240186816308773

# === Configurações Gerais de Qualidade ===
bg_color white
set ray_shadow, 0
set antialias, 2
set orthoscopic, 0
set valence, 0

# Centra e ajusta visualização
center 1TUP_sasa
zoom 1TUP_sasa
orient 1TUP_sasa

# Comandos úteis:
# hide surface - ocultar superfície
# hide sticks - ocultar bastões
# hide cartoon - ocultar cartoon
# set transparency, 0.7 - aumentar transparência

# Para salvar a sessão, use:
# save /Users/madsonluna/Downloads/biohub-alpha/1TUP_sasa.pse