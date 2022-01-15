const lightCodeTheme = require('prism-react-renderer/themes/github');
const darkCodeTheme = require('prism-react-renderer/themes/dracula');

// With JSDoc @type annotations, IDEs can provide config autocompletion
/** @type {import('@docusaurus/types').DocusaurusConfig} */
(module.exports = {
  title: 'fenics-optimize',
  tagline: 'fenics-optimize enables reusable and straightforward UFL coding for physical optimization problems and provides decorators that bridge easily between a fenics calculation chain and optimizers.',
  url: 'https://fenics-optimize.naruki-ichihara.com',
  baseUrl: '/',
  onBrokenLinks: 'throw',
  onBrokenMarkdownLinks: 'warn',
  favicon: 'img/favicon.ico',
  organizationName: 'Naruki-Ichihara', // Usually your GitHub org/user name.
  projectName: 'fenics_optimize_examples', // Usually your repo name.

  presets: [
    [
      '@docusaurus/preset-classic',
      /** @type {import('@docusaurus/preset-classic').Options} */
      ({
        docs: {
          sidebarPath: require.resolve('./sidebars.js'),
          // Please change this to your repo.
          editUrl: 'https://github.com/facebook/docusaurus/edit/main/website/',
        },
        blog: {
          showReadingTime: true,
          // Please change this to your repo.
          editUrl:
            'https://github.com/facebook/docusaurus/tree/main/packages/create-docusaurus/templates/shared/',
        },
        theme: {
          customCss: require.resolve('./src/css/custom.css'),
        },
      }),
    ],
  ],

  themeConfig:
    /** @type {import('@docusaurus/preset-classic').ThemeConfig} */
    ({
      navbar: {
        title: 'fenics-optimize',
        logo: {
          alt: 'My Site Logo',
          src: 'img/logo.svg',
        },
        items: [
          {
            type: 'doc',
            docId: 'intro',
            position: 'left',
            label: 'Documentation',
          },
          {
            type: 'doc',
            docId: 'intro',
            position: 'left',
            label: 'Blog',
          },
          {
            href: 'https://github.com/Naruki-Ichihara/fenics_optimize_examples',
            label: 'GitHub',
            position: 'right',
          },
        ],
      },
      footer: {
        style: 'light',
        links: [
          {
            title: 'Docs',
            items: [
              {
                label: 'Instllation',
                to: '/docs/intro',
              },
              {
                label: 'Getting started',
                to: '/docs/intro',
              },
              {
                label: 'Cookbook',
                to: '/docs/intro',
              }
            ],
          },
        ],
        copyright: `Copyright © ${new Date().getFullYear()} Naruki Ichihara. `,
      },
      prism: {
        theme: lightCodeTheme,
        darkTheme: darkCodeTheme,
      },
    }),
});