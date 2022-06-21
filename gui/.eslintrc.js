module.exports = {
    "extends": [
        "standard",
        "plugin:react/recommended"
    ],
    "parser": "babel-eslint",
    "parserOptions": {
        "ecmaFeatures": {
            "jsx": true,
            "modules": true
        }
    },
    "globals": {
        "fetch": false,
        "browser": true
    },
    "plugins": [
        "react", "react-hooks", "testing-library", "jest"
    ],
    "rules": {
        "space-before-function-paren": ["error", {
            "asyncArrow": "always",
            "named": "never",
            "anonymous": "never"
        }],
        "camelcase": [0],
        "no-async-promise-executor": ["warn"],
        "object-curly-spacing": [0],
        "object-curly-newline": [0],
        "no-var": [0],
        "prefer-const": [0],
        "quotes": [0],
        "quote-props": [0],
        "indent": [0],
        "multiline-ternary": [0],
        "no-empty": [0],
        "object-shorthand": [0],
        "dot-notation": [0],
        "lines-between-class-members": [0],
        "react-hooks/rules-of-hooks": "warn",
        "react-hooks/exhaustive-deps": "warn",
        "react/display-name": [0]
    },
    "settings": {
        "react": {
            "version": "detect"
        }
    },
    "overrides": [
        {
            "files": [
                "**/*.spec.js",
                "**/*.spec.jsx"
            ],
            "env": {
                "jest": true
            }
        }
    ]
}