// VENDORED from DocumenterVitepress (template/src/.vitepress/config.mts).
//
// DocumenterVitepress substitutes its own copy only when this file is absent, so
// this copy is frozen at the version it was taken from. There are two local
// changes, both marked below: the favicon links (DocumenterVitepress emits a
// single hardcoded favicon.ico link and offers no hook for a second, so it
// cannot follow the reader's light/dark browser theme), and the navbar logo.
//
// On a DocumenterVitepress upgrade: diff this against the new template and
// re-apply that one line. Every placeholder token in this file must stay intact
// -- the build substitutes them in place by exact string match.

import { defineConfig } from 'vitepress'
import { tabsMarkdownPlugin } from 'vitepress-plugin-tabs'
import { mathjaxPlugin } from './mathjax-plugin'
import { juliaReplTransformer } from './julia-repl-transformer'
import footnote from "markdown-it-footnote";
import path from 'path'

const mathjax = mathjaxPlugin()

function getBaseRepository(base: string): string {
  if (!base || base === '/') return '/';
  const parts = base.split('/').filter(Boolean);
  return parts.length > 0 ? `/${parts[0]}/` : '/';
}

const baseTemp = {
  base: 'REPLACE_ME_DOCUMENTER_VITEPRESS',// TODO: replace this in makedocs!
}

const navTemp = {
  nav: 'REPLACE_ME_DOCUMENTER_VITEPRESS',
}

const nav = [
  ...navTemp.nav,
  {
    component: 'VersionPicker'
  }
]

// https://vitepress.dev/reference/site-config
export default defineConfig({
  base: 'REPLACE_ME_DOCUMENTER_VITEPRESS',// TODO: replace this in makedocs!
  title: 'REPLACE_ME_DOCUMENTER_VITEPRESS',
  description: 'REPLACE_ME_DOCUMENTER_VITEPRESS',
  lastUpdated: true,
  cleanUrls: true,
  outDir: 'REPLACE_ME_DOCUMENTER_VITEPRESS', // This is required for MarkdownVitepress to work correctly...
  head: [
    // Dark-mode favicon, shipped as two .ico files rather than SVG. Safari ranks
    // an .ico above an SVG icon and paints the .ico even when both are offered,
    // so an SVG variant would never be seen there. Two further mechanisms that
    // look like they should work do not either: a `prefers-color-scheme` rule
    // inside the icon is evaluated as light whatever the OS theme, and `media`
    // on the link is not honoured everywhere. So the media-scoped links are the
    // progressive-enhancement base and the script is the fallback -- it removes
    // and re-appends rather than assigning href, because changing href alone
    // does not always make a browser re-read the icon, and it re-runs when the
    // OS theme changes. The DocumenterVitepress favicon placeholder is dropped:
    // it only ever emits one fixed `favicon.ico` link, which Safari would then
    // prefer over these.
    ['link', { rel: 'icon', media: '(prefers-color-scheme: light)', href: `${baseTemp.base}favicon-light.ico` }],
    ['link', { rel: 'icon', media: '(prefers-color-scheme: dark)',  href: `${baseTemp.base}favicon-dark.ico` }],
    ['script', {}, `(function(){var b=${JSON.stringify(baseTemp.base)},q=window.matchMedia('(prefers-color-scheme: dark)');function a(){var o=document.querySelectorAll('link[rel="icon"]');for(var i=0;i<o.length;i++)o[i].parentNode.removeChild(o[i]);var l=document.createElement('link');l.rel='icon';l.href=b+(q.matches?'favicon-dark.ico':'favicon-light.ico');document.head.appendChild(l);}q.addEventListener?q.addEventListener('change',a):q.addListener(a);a();})();`],
    ['script', {src: `${getBaseRepository(baseTemp.base)}versions.js`}],
    // ['script', {src: '/versions.js'], for custom domains, I guess if deploy_url is available.
    ['script', {src: `${baseTemp.base}siteinfo.js`}],
    // REPLACE_ME_DOCUMENTER_VITEPRESS_NOINDEX
  ],
  
  markdown: {
    codeTransformers: [juliaReplTransformer()],
    config(md) {
      md.use(tabsMarkdownPlugin);
      md.use(footnote);
      mathjax.markdownConfig(md);
    },
    theme: {
      light: "github-light",
      dark: "github-dark"
    },
  },
  vite: {
    plugins: [
      mathjax.vitePlugin,
    ],
    define: {
      __DEPLOY_ABSPATH__: JSON.stringify('REPLACE_ME_DOCUMENTER_VITEPRESS_DEPLOY_ABSPATH'),
    },
    resolve: {
      alias: {
        '@': path.resolve(__dirname, '../components')
      }
    },
    optimizeDeps: {
      exclude: [ 
        '@nolebase/vitepress-plugin-enhanced-readabilities/client',
        'vitepress',
        '@nolebase/ui',
      ], 
    }, 
    ssr: { 
      noExternal: [ 
        // If there are other packages that need to be processed by Vite, you can add them here.
        '@nolebase/vitepress-plugin-enhanced-readabilities',
        '@nolebase/ui',
      ], 
    },
  },
  themeConfig: {
    outline: 'deep',
    // Navbar logo, set directly instead of through the placeholder.
    // DocumenterVitepress only ever emits a single `logo: { src: '/logo.svg' }`,
    // but themeConfig.logo is a ThemeableImage, so VitePress can pick between the
    // two variants with its own `.dark` class -- no second asset that adapts via
    // prefers-color-scheme, and no CSS repaint. VPImage applies withBase(), so the
    // absolute paths survive a deploy under a base prefix. The build logs a
    // "No logo.png file found" warning because no logo.svg exists any more; it is
    // expected, and the replacement it would have made is a no-op here.
    logo: { light: '/logo-light.svg', dark: '/logo-dark.svg', width: 24, height: 24 },
    search: {
      provider: 'local',
      options: {
        detailedView: true
      }
    },
    nav,
    sidebar: 'REPLACE_ME_DOCUMENTER_VITEPRESS',
    sidebarDrawer: 'REPLACE_ME_DOCUMENTER_VITEPRESS_SIDEBAR_DRAWER',
    editLink: 'REPLACE_ME_DOCUMENTER_VITEPRESS',
    socialLinks: [
      { icon: 'github', link: 'REPLACE_ME_DOCUMENTER_VITEPRESS' }
    ],
    footer: {
      message: 'Made with <a href="https://luxdl.github.io/DocumenterVitepress.jl/dev/" target="_blank"><strong>DocumenterVitepress.jl</strong></a><br>',
      copyright: `© Copyright ${new Date().getUTCFullYear()}.`
    }
  }
})
