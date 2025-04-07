// tailwind.config.js
/** @type {import('tailwindcss').Config} */
module.exports = {
    content: [
      "./src/**/*.{js,jsx,ts,tsx}",
      "./public/index.html",
    ],
    darkMode: 'class',
    theme: {
      extend: {
        colors: {
          gray: {
            750: '#2b3544',
            850: '#1a2231',
            950: '#0f1521',
          },
          blue: {
            350: '#7dabf8',
            450: '#4f85e6',
            550: '#3a72d6',
            650: '#2a5db8',
            750: '#214a98',
            850: '#183878',
          },
        },
        spacing: {
          '72': '18rem',
          '84': '21rem',
          '96': '24rem',
          '128': '32rem',
        },
        maxWidth: {
          '8xl': '90rem',
        },
        opacity: {
          '15': '0.15',
          '35': '0.35',
          '85': '0.85',
          '95': '0.95',
        },
        borderRadius: {
          'xl': '1rem',
          '2xl': '2rem',
          '3xl': '3rem',
        },
        boxShadow: {
          'inner-lg': 'inset 0 2px 4px 0 rgba(0, 0, 0, 0.2)',
          'inner-xl': 'inset 0 4px 8px 0 rgba(0, 0, 0, 0.25)',
          'soft-xl': '0 10px 25px -5px rgba(0, 0, 0, 0.1), 0 10px 10px -5px rgba(0, 0, 0, 0.04)',
          'blue-glow': '0 0 15px rgba(66, 153, 225, 0.5)',
        },
        height: {
          '128': '32rem',
          '144': '36rem',
        },
        minHeight: {
          '24': '6rem',
          '32': '8rem',
          '48': '12rem',
          '64': '16rem',
          '96': '24rem',
        },
        fontFamily: {
          sans: [
            'Inter',
            'ui-sans-serif',
            'system-ui',
            '-apple-system',
            'BlinkMacSystemFont',
            '"Segoe UI"',
            'Roboto',
            '"Helvetica Neue"',
            'Arial',
            '"Noto Sans"',
            'sans-serif',
          ],
          mono: [
            'JetBrains Mono',
            'ui-monospace',
            'SFMono-Regular',
            'Menlo',
            'Monaco',
            'Consolas',
            '"Liberation Mono"',
            '"Courier New"',
            'monospace',
          ],
        },
        animation: {
          'pulse-slow': 'pulse 3s cubic-bezier(0.4, 0, 0.6, 1) infinite',
          'pulse-fast': 'pulse 1s cubic-bezier(0.4, 0, 0.6, 1) infinite',
          'bounce-slow': 'bounce 2s infinite',
        },
        zIndex: {
          '60': '60',
          '70': '70',
          '80': '80',
          '90': '90',
          '100': '100',
        },
        backdropBlur: {
          'xs': '2px',
        },
        typography: {
          DEFAULT: {
            css: {
              color: '#e2e8f0',
              a: {
                color: '#3b82f6',
                '&:hover': {
                  color: '#60a5fa',
                },
              },
              h1: {
                color: '#f8fafc',
              },
              h2: {
                color: '#f8fafc',
              },
              h3: {
                color: '#f8fafc',
              },
              h4: {
                color: '#f8fafc',
              },
              strong: {
                color: '#f1f5f9',
              },
              code: {
                color: '#f1f5f9',
                backgroundColor: '#1e293b',
                borderRadius: '0.25rem',
                padding: '0.2em 0.4em',
              },
              pre: {
                backgroundColor: '#0f172a',
                color: '#f8fafc',
              },
              blockquote: {
                color: '#cbd5e1',
                borderLeftColor: '#475569',
              },
            },
          },
        },
      },
    },
    plugins: [
      require('@tailwindcss/forms'),
      require('@tailwindcss/typography'),
    ],
  };