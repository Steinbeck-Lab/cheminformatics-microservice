import { defineConfig } from 'vitepress'

// https://vitepress.dev/reference/site-config
export default defineConfig({
  title: "Cheminformatics Microservice",
  description: "This set of essential and valuable microservices is designed to be accessed via API calls to support cheminformatics. Generally, it is designed to work with SMILES-based inputs and could be used to translate between different machine-readable representations, get Natural Product (NP) likeliness scores, visualize chemical structures, and generate descriptors. In addition, the microservice also host an instance of STOUT and another instance of DECIMER (two deep learning models for IUPAC name generation and optical chemical structure recognition, respectively).",

  themeConfig: {
    search: {
      provider: 'local'
    },
    
    logo: { 
      light: 'logo.png',
      dark: 'logo_light.png'
    },

    siteTitle: '',
    
    nav: [
      { text: 'Home', link: '/' },
      { text: 'Docs', link: '/introduction' },
      { text: 'API', link: 'https://api.naturalproducts.net/latest/docs' },
    ],

    sidebar: [
      {
        text: 'Welcome',
        items: [
          { text: 'Introduction', link: '/introduction' },
          { text: 'Versions', link: '/versions' },
          { text: 'Architecture', link: '/architecture' },
          { text: 'Sustainability', link: '/sustainability' },
        ]
      },
      {
        text: 'Installation',
        items: [
          { text: 'Docker', link: '/docker' },
          { text: 'Cluster Deployment (K8S)', link: '/cluster-deployment' },
          { text: 'Standalone (Python virtual environment)', link: '/standalone' },
          { text: 'Scaling', link: '/scaling' },
        ]
      },
      {
        text: 'Modules',
        items: [
          { text: 'Chem', link: '/chem' },
          { text: 'Convert', link: '/conversions' },
          { text: 'Depict', link: '/depict' },
          { text: 'Tools', link: '/tools' },
          { text: 'OCSR', link: '/ocsr' },
        ]
      },
      {
        text: 'Development',
        items: [
          { text: 'Local Setup', link: '/installation' },
          { text: 'License', link: '/license' },
          { text: 'Issues/Feature requests', link: '/issues' },
          { text: 'Contributors', link: '/contributors' }
        ]
      }
    ],

    socialLinks: [
      { icon: 'github', link: 'https://github.com/Steinbeck-Lab/cheminformatics-microservice' }
    ],

    footer: {
      message: 'Source code released under the MIT License | Data are provided under the Creative Commons Attribution (aka CC-BY 4.0) <br/> Funded by the <a href="https://www.dfg.de/en/index.jsp" style="color: blue" target="_blank">Deutsche Forschungsgemeinschaft (DFG, German Research Foundation)</a> under the <a href="https://www.nfdi4chem.de/" style="color: blue" target="_blank">National Research Data Infrastructure – NFDI4Chem</a> – Project number: 441958208 and <a href="https://www.chembiosys.de/en/research/projects/project-inf.html" style="color: green" target="_blank">ChemBioSys (Project INF)</a> - Project number: 239748522 - SFB 1127.',
      copyright: `Copyright © ${new Date().getFullYear()} Steinbeck Lab`
    }
  }
})
