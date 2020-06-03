module.exports = {
    "extends": [
        "standard",
        "plugin:react/recommended",
        "plugin:react-hooks/recommended"
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
    "rules": {
        "space-before-function-paren": ["error", "never"],
        "camelcase": [0],
        "react-hooks/exhaustive-deps": "off"
    },
    "settings": {
        "react": {
            "version": "detect"
        }
    }
}