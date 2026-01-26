#quarto publish quarto-pub --no-render  
# in r run
rsconnect::deployApp(appDir = '~/analysis/flow/w0_analysis/w0_book/_book',recordDir = '~/analysis/flow/w0_analysis/w0_book', account = 'mgarofal')
