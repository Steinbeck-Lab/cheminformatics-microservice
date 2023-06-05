import { defineConfig } from 'vitepress'

// https://vitepress.dev/reference/site-config
export default defineConfig({
  title: "Cheminformatics Python Microservice",
  description: "This set of essential and valuable microservices is designed to be accessed via API calls to support cheminformatics. Generally, it is designed to work with SMILES-based inputs and could be used to translate between different machine-readable representations, get Natural Product (NP) likeliness scores, visualize chemical structures, and generate descriptors. In addition, the microservices also host an instance of STOUT and another instance of DECIMER (two deep learning models for IUPAC name generation and optical chemical structure recognition, respectively).",
  base: '/cheminformatics-python-microservice/',
  themeConfig: {
    logo: { 
      light: 'logo.png',
      dark: 'logo_light.png'
    },

    siteTitle: '',
    
    // https://vitepress.dev/reference/default-theme-config
    nav: [
      { text: 'Home', link: '/' },
      { text: 'Docs', link: '/introduction' },
      { text: 'API', link: 'https://api.naturalproducts.net/docs' }
    ],

    sidebar: [
      {
        text: 'Welcome',
        items: [
          { text: 'Introduction', link: '/introduction' },
          { text: 'Versions', link: '/versions' },
          { text: 'Architecture', link: '/architecture' },
        ]
      },
      {
        text: 'Usage',
        items: [
          { text: 'Public API', link: '/public-api' },
          { text: 'Docker', link: '/docker' },
          { text: 'Cluster Deployment (K8S)', link: '/cluster-deployment' }
        ]
      },
      {
        text: 'Modules',
        items: [
          { text: 'Conversions', link: '/conversions' },
          { text: 'Chem', link: '/chem' },
          { text: 'DECIMER', link: '/decimer' },
          { text: 'STOUT', link: '/stout' }
        ]
      },
      {
        text: 'Development',
        items: [
          { text: 'Local Installation', link: '/installation' },
          { text: 'License', link: '/license' },
          { text: 'Issues/Feature requests', link: '/issues' },
          { text: 'Contributors', link: '/contributors' }
        ]
      }
    ],

    socialLinks: [
      { icon: 'github', link: 'https://github.com/Steinbeck-Lab/cheminformatics-python-microservice' }
    ]
  }
})
