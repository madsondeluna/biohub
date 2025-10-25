# Script PyMOL gerado pelo BioHub
# Propriedade: sasa

# Carrega a estrutura
load /Users/madsonluna/Documents/biohub-1/4HHB_sasa.pdb, 4HHB_sasa

# Remove todas as representações padrão
hide everything, 4HHB_sasa

# === SASA (Acessibilidade ao Solvente) ===
# Vermelho = Enterrado (0.0 Ų), Azul = Exposto (10.8 Ų)

# Representação Cartoon (fita)
show cartoon, 4HHB_sasa
cartoon automatic, 4HHB_sasa
set cartoon_fancy_helices, 1
spectrum b, red_white_blue, 4HHB_sasa, minimum=0.0, maximum=10.788427721083055

# Representação Sticks (bastões)
show sticks, 4HHB_sasa
set stick_radius, 0.2, 4HHB_sasa
set stick_color, gray, 4HHB_sasa
util.cbag 4HHB_sasa

# Representação Surface (superfície)
show surface, 4HHB_sasa
set surface_quality, 1
set transparency, 0.5, 4HHB_sasa
# Aplica gradiente de SASA na superfície
set surface_color, white, 4HHB_sasa
spectrum b, red_white_blue, 4HHB_sasa, minimum=0.0, maximum=10.788427721083055

# === Configurações Gerais de Qualidade ===
bg_color white
set ray_shadow, 0
set antialias, 2
set orthoscopic, 0
set valence, 0

# Centra e ajusta visualização
center 4HHB_sasa
zoom 4HHB_sasa
orient 4HHB_sasa

# Comandos úteis:
# hide surface - ocultar superfície
# hide sticks - ocultar bastões
# hide cartoon - ocultar cartoon
# set transparency, 0.7 - aumentar transparência

# Para salvar a sessão, use:
# save /Users/madsonluna/Documents/biohub-1/4HHB_sasa.pse