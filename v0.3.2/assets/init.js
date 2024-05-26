require.config({
    paths: {
        mermaid: "https://cdn.jsdelivr.net/npm/mermaid/dist/mermaid"
    }
});
require(['mermaid'], function(mermaid) { mermaid.init() });